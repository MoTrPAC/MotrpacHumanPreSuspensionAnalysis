#'  @title SCION implementation for preCAWG.
#' Combines preprocessing and network inference into one function
#'
#' This function is based on the source code \code{https://github.com/nmclark2/SCION}
#'
#' @author Natalie M Clark, Christopher Jin
#'
#' @param regulators the set of hypothesized tissue/omes that regulate the
#'    targets matrix, from .load_scion_matrixes
#' @param targets the set of hypothesized tissue/omes that are targetted by the
#'    regulators matrix, from \code{.load_scion_matrixes}
#' @param randomGroupCode - regulators, targets, and clustering is all done on a
#'    per-modality basis. All the clustering is done using z-scores from the differential analysis
#' @param permute number of random permutations to perform for edge trimming.
#'    Default NULL, which means no permutations will be performed
#' @param dim dimension on which to permute. options are "row" or "col". default "col"
#' @param cluster boolean option to cluster using c-means prior to network inference. default TRUE
#' @param dir.name name of directory to save results. default "exp".
#'    if using permutations, this parameter is ignored, and the directories are
#'    named after the permutation number (1,2,3...)
#' @param weightthreshold threshold for edge trimming. all edges with weight < weightthreshold
#'    are removed from the network. minimum value of 0. default 0.
#' @param normalize boolean to normalize edge weights to a 0-1 scale. default FALSE
#' @param num.cores number of cores to use for parallelization. num.cores-1 will
#'     be used for parallelization. if num.cores < 3, no parallelization is performed. default 1.
#' @param connect.hubs boolean to connect the hubs between clusters.
#'    this parameter is ignored if clustering is not performed. default TRUE
#' @param verbose boolean to display detailed output. default FALSE
#'
#' @export run_SCION
#'
#' @importFrom stringr str_extract
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr mutate filter select

run_SCION <- function (randomGroupCode = c("ADUResist", "ADUEndur"),
                       regulators,
                       targets,
                       permute = NULL,
                       dim = "col",
                       cluster = TRUE,
                       dir.name = "exp",
                       weightthreshold = 0,
                       normalize = FALSE,
                       num.cores = 1,
                       connect.hubs = TRUE,
                       verbose = TRUE) {
  check_package_installation(pkg = "pacman")
  check_package_installation(pkg = "parallel")
  check_package_installation(pkg = "doParallel")
  check_package_installation(pkg = "randomForest")
  check_package_installation(pkg = "MotrpacHumanPreSuspensionData")

  message("Note: epigenetic targets are not currently supported")

  #create directory to save results
  if(!is.null(permute)){
    my.dir <- file.path(dir.name, as.character(permute))
    cat(paste("Starting permutation ",permute,"\n",sep=""))
  }else{
    my.dir <- dir.name
  }
  message(my.dir)
  if(!dir.exists(my.dir)) dir.create(my.dir, recursive = TRUE)

  cat("Processing data tables\n") #first need to pre-process the data
  tables = scion_data_processing(randomGroupCode = randomGroupCode,
                                 reg.mat = regulators,
                                 target.mat = targets,
                                 permute,
                                 dim)
  #Here instead of directly entering in the prefix at function start, I just
  #assign it according to the regulator syntax.
  #I just choose the first ome that shows up in the regulators in the first row
  prefix = unlist(strsplit(rownames(tables[["reg"]])[1], "..", fixed = TRUE))[2]
  if(grepl("metab", prefix)) prefix = "metab"

  #assign pre-processed data tables
  mytargetdata <- as.data.frame(tables$target)
  myregdata <- as.data.frame(tables$reg)
  #make.names on rownames to ensure no future errors
  rownames(mytargetdata) <- make.names(rownames(mytargetdata))
  rownames(myregdata) <- make.names(rownames(myregdata))

  #We need to manipulate the input structure to run properly.
  if(cluster){
    reg_tissue_omes = lapply(strsplit(rownames(myregdata), "\\.{3}|(?<!\\.)\\.\\.", perl=TRUE),
                             function(x) if (length(x) > 3) tail(x, 3) else x)

    target_tissue_omes = strsplit(rownames(mytargetdata), "\\.\\.")

    split_rownames_df = as.data.frame(do.call(rbind, c(reg_tissue_omes, target_tissue_omes)))
    colnames(split_rownames_df) = c("feature_id", "ome", "tissue")
    unique_pairs = split_rownames_df %>%
      dplyr::mutate(ome = dplyr::case_when(
        grepl("metab", ome) ~ "metab",
        TRUE ~ ome
      )) %>%
      dplyr::distinct(ome, tissue)
    #----building expected DA_list format for run_cmeans--------
    DA_Input_Cmeans = list()
    for(row in seq(nrow(unique_pairs))){
      selected_ome = unique_pairs[row,][["ome"]]
      selected_tissue = unique_pairs[row,][["tissue"]]
      selected_ome_mod = gsub("\\.", "-", selected_ome) #have to do this bc the rownames arent happy with "-"
      if(selected_tissue == "adipose" && (selected_ome_mod == "prot-pr" || selected_ome_mod == "prot-ph"))
        stop("Adipose Prot-ph, Prot-pr only have 2 timepoints. We don't recommend using these, it messes with the clustering")
      da_table = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(selected_omes = selected_ome_mod,
                                                                           selected_tissues = selected_tissue,
                                                                           single_matrix = TRUE)
      if(randomGroupCode == "ADUEndur") da_table = da_table %>%
        dplyr::filter(contrast_type == "exercise_with_controls") %>%
        dplyr::filter(contrast_category == "EE-CON")
      if(randomGroupCode == "ADUResist") da_table = da_table %>%
        dplyr::filter(contrast_type == "exercise_with_controls") %>%
        dplyr::filter(contrast_category == "RE-CON")

      filt_to_inputs = split_rownames_df %>%
        dplyr::mutate(ome = dplyr::case_when(
          grepl("metab", ome) ~ "metab",
          TRUE ~ ome
        )) %>%
        dplyr::filter(ome == selected_ome, tissue == selected_tissue) %>%
        dplyr::pull(feature_id)
      #i add "__" to make splitting the results and manipulating syntax a bit easier
      da_filt = da_table %>%
        dplyr::filter(feature_id %in% filt_to_inputs) %>%
        dplyr::mutate(feature_id = paste(feature_id, selected_tissue, sep = "__"))

      #expected format for c-means, cameraPR unnesting
      merged_name = paste(selected_tissue, selected_ome_mod, sep = ".")
      DA_Input_Cmeans[[merged_name]] = da_filt
    }
    clusters = scion_run_cmeans(DA_Input_Cmeans)

    cmeans_results = clusters[["cluster"]] #again we force everything to be one given tissue
    #now we have to manipulate the names of these back into the same format as mytarget/reg data
    cluster_org = data.frame(name = names(cmeans_results), cluster = as.integer(cmeans_results)) %>%
      dplyr::mutate(name = gsub("-", ".", name)) %>%
      dplyr::mutate(
        ome = str_extract(name, "^[^ ]+"),                        # Value before the space
        tissue = str_extract(name, "(?<=__)[^_]+$"),              # Value after the last "__"
        feature_id = str_extract(name, "(?<= )[A-Za-z0-9._]+(?=__)") # Value after space and before "__"
      ) %>%
      dplyr::mutate(rowname_format = paste(feature_id, ome, tissue, sep = "..")) %>%
      dplyr::select(rowname_format, cluster) %>%
      tibble::column_to_rownames("rowname_format")
  }

  #infer a network using GENIE3 on each cluster
  SCION_infer <- function(mytargetdata,
                          myregdata,
                          clusterresults,
                          weightthreshold,
                          prefix,
                          permute=NULL,
                          normalize=FALSE,
                          verbose=F){

    finalnetwork <- data.frame(Regulator=character(),
                               Interaction=character(),
                               Target=character(),
                               Weight=double(),
                               Cluster=numeric(),
                               stringsAsFactors=FALSE)
    myhubs=NULL
    cat("Inferring cluster-specific networks.\n")
    for (i in 1:max(clusterresults$cluster)){
      #so here we need to reverse back into selecting the tissue
      mygenes = row.names(clusterresults)[clusterresults$cluster==i]
      clustertargetdata = mytargetdata[row.names(mytargetdata)%in%mygenes,]
      clusterregdata = myregdata[row.names(myregdata)%in%mygenes,]

      #we need at least one target and at least one regulator
      if (dim(clustertargetdata)[1]<1 | dim(clusterregdata)[1]<1){
        next
      }

      #infer the network
      if(verbose){
        cat(paste("Inferring network for cluster ",i,
                  " with ", dim(clusterregdata)[1], " regulators and ",
                  dim(clustertargetdata)[1], " targets.\n",sep=""))
      }
      network = RS.Get.Weight.Matrix(t(clustertargetdata),t(clusterregdata),normalize=normalize,num.cores=num.cores)
      #if network inference failed, move on
      if (is.null(network)){
        next
      }

      #now make a new network where we eliminate all the low confidence edges
      trimmednet = data.frame(ifelse(network<weightthreshold,NaN,network))

      #translate the trimmed network into a table we can import into cytoscape
      networktable = data.frame(Regulator=character(), Interaction=character(), Target=character(), Weight=double(),Cluster=numeric(),
                                stringsAsFactors=FALSE)
      row=1
      for (j in 1:dim(trimmednet)[1]){
        for (k in 1:dim(trimmednet)[2]){
          #skip NaNs as these have no edge
          if (is.na(trimmednet[j,k])){
            next
          }else{
            networktable[row,] = cbind(colnames(trimmednet)[k],"regulates",rownames(trimmednet)[j],trimmednet[j,k],i)
            row = row+1
          }
        }
      }

      #save this network
      finalnetwork = rbind(finalnetwork,networktable)
      write.csv(finalnetwork,paste0(my.dir,"/network_cluster_",i,".csv"))

      #get the hub gene (highest outdegree) and save it to connect the clusters later
      #if there is a tie, we save both hubs
      myregs = unique(networktable$Regulator)
      for (j in 1:length(myregs)){
        numedges = sum(networktable$Regulator%in%myregs[j])
        if (j==1){
          hubedges = numedges
          hub = myregs[j]
        }else if (numedges>hubedges){
          hubedges = numedges
          hub = myregs[j]
        }else if (numedges==hubedges){
          hub = rbind(hub,myregs[j])
        }
      }
      if (!exists("myhubs")){
        myhubs = hub
      } else{
        myhubs = rbind(myhubs,hub)
      }
    }

    #connect the hubs for each cluster
    if (connect.hubs && exists("myhubs") && length(myhubs)>2){
      cat("Connecting hubs\n")
      write.csv(myhubs,paste0(my.dir,"/hubs.csv"))
      write.csv(mytargetdata,paste0(my.dir,"/mytargetdata.csv"))
      write.csv(myregdata,paste0(my.dir,"/myregdata.csv"))

      hubtargetdata = mytargetdata[row.names(mytargetdata)%in%myhubs,]
      hubregdata = myregdata[row.names(myregdata)%in%myhubs,]

      #Note: this pulls independent of tissue any features that match up to
      #gene symbols related to prot-ph
      if (dim(hubtargetdata)[1]==0 & prefix=='prot-ph'){
        #strip the PTM information so that we can get the targets
        genes <- unlist(strsplit(myhubs,'..', fixed = TRUE))
        genes <- genes[seq(1,length(genes),by=3)]

        any_ome_ids = MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
          dplyr::filter(feature_id %in% genes) %>%
          dplyr::left_join(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE,
                           by = "gene_symbol") %>%
          dplyr::filter(!assay.y %in% c("epigen-methylcap-seq", "epigen-atac-seq")) %>%
          dplyr::pull(feature_id.y)
        #----
        hubtargetdata = mytargetdata %>%
          dplyr::filter(sub("\\.\\..*", "", rownames(.)) %in% any_ome_ids)

      }
      else if (dim(hubtargetdata)[1]==0 & prefix=="metab"){
        #use the regulator matrix AS the target matrix for these
        #good for metabolome
        hubtargetdata = hubregdata
      }
      network = RS.Get.Weight.Matrix(t(hubtargetdata),t(hubregdata),normalize=normalize)
      if (!is.null(network)){
        #now make a new network where we eliminate all the low confidence edges
        #There are some network values with <0 that do get trimmed
        trimmednet = data.frame(ifelse(network<weightthreshold,NaN,network))
        #translate the trimmed network into a table we can import into cytoscape
        networktable = data.frame(Regulator=character(), Interaction=character(), Target=character(), Weight=double(),Cluster=numeric(),
                                  stringsAsFactors=FALSE)
        row=1
        for (j in 1:dim(trimmednet)[1]){
          for (k in 1:dim(trimmednet)[2]){
            #skip NaNs as these have no edge
            if (is.na(trimmednet[j,k])){
              next
            }else{
              networktable[row,] = cbind(colnames(trimmednet)[k],"regulates",rownames(trimmednet)[j],trimmednet[j,k],"")
              row = row+1
            }
          }
        }
        #save this network
        finalnetwork = rbind(finalnetwork,networktable)
      }
    }

    #save the final network and the clustering information
    if(!is.null(permute)){
      cat(paste("Saving network for permutation ",permute,"\n",sep=""))
      write.table(finalnetwork,paste (my.dir,"/",prefix, '-SCION-network-permutation-',permute,'.tsv', sep=''),row.names=FALSE,quote=FALSE,sep='\t')
    }else{
      write.table(finalnetwork,paste (my.dir, "/", prefix, '-SCION-network-full.tsv', sep=''),row.names=FALSE,quote=FALSE,sep='\t')
    }
  }

  #infer the network
  #set seed so network inference is the same each time
  set.seed(0)
  SCION_infer(mytargetdata,
              myregdata,
              clusterresults = cluster_org,
              weightthreshold,
              prefix,
              permute,
              normalize,
              verbose)
  cat("Done\n")
}


scion_data_processing <- function(randomGroupCode,
                                  reg.mat,
                                  target.mat,
                                  permute=NULL,
                                  dim="row"){
  if(is.character(reg.mat) | is.character(target.mat))
    stop("It looks like your input is a character. Use .load_scion_matrixes to specify your input more specifically")

  find_shared_participants = .find_shared_scion_matrixes(randomGroupCode, reg.mat, target.mat)
  reg.mat = find_shared_participants[["scion_matrix_one"]]
  target.mat = find_shared_participants[["scion_matrix_two"]]

  #need to change the seed for each permutation. otherwise, the shuffling is exactly the same when running in parallel
  if(!is.null(permute)){
    set.seed(permute)
  }
  #permute data matrices if desired
  #target mat
  if(!is.null(permute)){
    if(dim=="col"){
      #permute the target matrix column-wise (by sample)
      target.names <- row.names(target.mat)
      target.mat <- apply(target.mat,2,sample)
      row.names(target.mat) <- target.names
    }else{
      #permute the target matrix row-wise (by gene)
      target.mat <- t(apply(target.mat,1,sample))
    }
  }
  #reg mat
  if(!is.null(permute)){
    if(dim=="col"){
      #permute the regulator matrix col-wise (by sample)
      reg.names <- row.names(reg.mat)
      reg.mat <- apply(reg.mat,2,sample)
      row.names(reg.mat) <- reg.names
    }else{
      #permute the regulator matrix row-wise (by site)
      reg.mat <- t(apply(reg.mat,1,sample))
    }
  }
  return(list("reg"=reg.mat,
              "target"=target.mat))
}

# target.matrix and input.matrix (TFs) have samples as rows and genes as columns
RS.Get.Weight.Matrix<- function(target.matrix,
                                input.matrix,
                                K="sqrt",
                                nb.trees=10000,
                                importance.measure="%IncMSE",
                                seed=NULL,
                                trace=TRUE,
                                normalize=TRUE,
                                num.cores=1, ...){
  #require(parallel)
  # set random number generator seed if seed is given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # to be nice, report when parameter importance.measure is not correctly spelled
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
  }

  # normalize expression matrix
  target.matrix <- apply(target.matrix, 2, function(x) { (x - mean(x,na.rm = T)) / sd(x,na.rm = T) } )
  input.matrix <- apply(input.matrix, 2, function(x) { (x - mean(x,na.rm =T)) / sd(x,na.rm=T) } )
  input.matrix <- input.matrix[,!is.na(colSums(input.matrix))]
  # setup weight matrix
  num.samples <- dim(target.matrix)[1]
  num.targets <- dim(target.matrix)[2]
  num.inputs <- dim(input.matrix)[2]
  target.names <- colnames(target.matrix)
  input.names <- colnames(input.matrix)
  #print(input.names)
  #if no inputs or targets, return NULL
  if (is.null(num.inputs) | is.null(num.targets)){
    return(NULL)
  }

  weight.matrix <- matrix(0.0, nrow=num.targets, ncol=num.inputs)
  rownames(weight.matrix) <- target.names
  colnames(weight.matrix) <- input.names

  # set mtry
  if (is.numeric(K)) {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.inputs))
  } else if (K == "all") {
    mtry <- num.inputs-1
  } else {
    stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
  }
  # if (trace) {
  #   cat(paste("Starting RF computations with ", nb.trees,
  #             " trees/target gene,\nand ", mtry,
  #             " candidate input genes/tree node\n",
  #             sep=""))
  #   flush.console()
  # }

  # compute importances for every target gene
  names(target.names)<-target.names
  #return(list(target.names,input.matrix,target.matrix))

  #parallelize if at least 3 cores, otherwise, don't
  if(num.cores>2){
    clst <- parallel::makeCluster(num.cores-1,type="FORK",outfile="log.txt")
    registerDoParallel(clst)
    imList<-parLapply(cl=clst, X=target.names, function(x) RSGWM2(x,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...))
    stopCluster(cl=clst)
  }else{
    imList<-lapply(target.names,function(x) RSGWM2(x,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...))
  }

  #return(imList)
  for(nm in names(imList))
  {
    tcols<-names(imList[[nm]])
    weight.matrix[nm,tcols] <- imList[[nm]]
  }
  # weight.matrix<-sapply(imList,function(x) x)

  mynet <- weight.matrix/num.samples
  if(normalize==TRUE){
    mynet <- (mynet-min(mynet,na.rm=TRUE))/(max(mynet,na.rm=TRUE)-min(mynet,na.rm=TRUE))
  }
  return(mynet)
  #    return(list(Weight=weight.matrix,Model=model.matrix,PedictionCorrelations=cor.vec))
}

RSGWM2<-function(target.gene.name,
                 num.targets,
                 target.names,
                 input.matrix,
                 target.matrix,
                 trace,
                 mtry,
                 nb.trees,
                 importance.measure,
                 ...){
  target.gene.idx<-which(target.names==target.gene.name)
  if (trace)
  {
    #cat(paste("Computing gene ", target.gene.idx, "/", num.targets, "\n", sep=""))
    flush.console()
  }
  #target.gene.name <- target.names[target.gene.idx]
  # remove target gene from input genes
  #these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
  temp.input.matrix<-input.matrix

  #do NOT remove the gene from the input matrix
  #this breaks the code when there is only 1 regulator in a network
  #I believe this was originally written to remove autoregulation from the network

  # if(target.gene.name %in% colnames(input.matrix))
  # {
  #   rmind<-which(colnames(input.matrix)==target.gene.name)
  #   #rmind<-grep(target.gene.name,colnames(input.matrix),fixed=T)
  #   temp.input.matrix<-input.matrix[,-rmind]
  #   #print(paste("Removing:",target.gene.name,"fom Input | Index Number:",rmind," |  New Dimensions:",dim(temp.input.matrix)[1],"X",dim(temp.input.matrix)[2]))
  # }
  x <- temp.input.matrix
  y <- target.matrix[,target.gene.name]
  #incSamp<-names(y)[which(!is.na(y))]
  #x<-x[incSamp,]
  #y<-y[incSamp]

  rf <- randomForest(x = x, y = y, mtry=mtry, ntree=nb.trees, keep.forest=F, importance=TRUE,...)

  im <- importance(rf)[,importance.measure]

  #im.names <- names(im)
  return(im)
}

#' @title Load data input for scion
#' @description Function to load data that fits the requirements of the input
#'    regulators and targets
#'
#' @param desired_matrixes hypothesized regulators/targets, which has to be in the order of tissue.ome using available values from
#' `tissue_available_list()` and `ome_available_list`, seperated by a "." (e.g. "blood.transcript-rna-seq")
#' You can use "metabolomics", which will include all metabolomics platforms
#' @return a list with targets and regulator matrixes
#'
#' @export .load_scion_matrixes
#'
#' @noRd

.load_scion_matrixes = function(desired_matrixes,
                                subset_TFs = TRUE,
                                subset_DE = TRUE){
  final_combined_matrix = data.frame()
  for(desired_input in desired_matrixes){
    desired_tissue = sub("\\..*", "", desired_input)
    if(!desired_tissue %in% tissue_available_list())
      stop("The syntax for regulators or targets isn't right. Regulators and targets should be in tissue.ome format (e.g blood.transcript-rna-seq)")
    desired_ome = sub(".*?\\.", "", desired_input)
    if(desired_ome == "metabolomics" | desired_ome == "metab") desired_ome = metab_only_list()
    if(!all(desired_ome %in% ome_available_list()))
      stop("The syntax for regulators or targets isn't right. Regulators and targets should be in tissue.ome format (e.g blood.transcript-rna-seq)")

    desired_qc_norm_single = MotrpacHumanPreSuspensionData::load_qc(selected_omes = desired_ome,
                                                                    selected_tissues = desired_tissue,
                                                                    remove_unnamed_metab = TRUE)
    #-------------------------------------------------------------------------------
    # here we handle parameters and use imputed matrixes for prot-(pr/ph)
    if(any(desired_ome == "prot-ph" | desired_ome == "prot-pr")) #have to add 'any' to handle multiple metab platforms
      desired_qc_norm_single[[desired_tissue]][[desired_ome]][["qc_norm"]] = desired_qc_norm_single[[desired_tissue]][[desired_ome]][["qc_imputed"]]
    #basically just copy over imputed into the name "qc_norm" because then I can just use the combine_qc_matrixes function as usual
    if(any(desired_ome == "prot-ph") & subset_TFs) desired_qc_norm_single = subset_qc(desired_qc_norm_single,
                                                                                      desired_features = MotrpacHumanPreSuspensionAnalysis::UTORONTO_TFs$feature_id)

    if(subset_DE){
      DE_features = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(selected_omes = desired_ome,
                                                                              selected_tissues = desired_tissue,
                                                                              single_matrix = TRUE) %>%
        dplyr::filter(contrast_type == "exercise_with_controls",
                      adj_p_value < 0.05) %>%
        dplyr::pull(feature_id)
      desired_qc_norm_single = subset_qc(desired_qc_norm_single, desired_features = DE_features)
    }
    #--------------------------------------------------------------------------------
    participant_format = combine_qc_matrixes(desired_qc_norm_single,
                                             make_metab_rownames = TRUE)
    #just make all the metab stuff just = "metab"

    if(nrow(final_combined_matrix) == 0){
      final_combined_matrix = participant_format
    }else{
      shared_cols = intersect(colnames(final_combined_matrix), colnames(participant_format))
      if(length(shared_cols) == 0 || any(is.na(shared_cols))) stop("desired targets or regulators have no shared participants")
      final_combined_matrix = final_combined_matrix %>% dplyr::select(dplyr::all_of(shared_cols))
      participant_format = participant_format %>% dplyr::select(dplyr::all_of(shared_cols))
      final_combined_matrix = rbind(final_combined_matrix, participant_format)
    }
  }
  unique_parts = sub("\\.\\..*", "", colnames(final_combined_matrix))
  message(paste("Matrix of", paste(desired_matrixes, collapse = ", "),
                "has", length(unique(unique_parts)), "participants,", length(unique_parts), "samples"))
  return(final_combined_matrix)
}


#' @title Find shared participant/timepoints for scion matrixes
#'
#' @details Because different omes and tissues have different participants included
#'    in some cases, this will filter and find shared columns between the two.
#'    Very similar to the implementation from \code{combine_qc_matrixes()}
#'
#' @param scion_matrix_one First matrix from \code{.load_scion_matrixes}
#' @param scion_matrix_two Second matrix from \code{.load_scion_matrixes}
#'
#' @return Both matrixes except subset to the same matching participants/timepoints
#' @export .find_shared_scion_matrixes
#'
#' @noRd
.find_shared_scion_matrixes = function(randomGroupCode,
                                       scion_matrix_one,
                                       scion_matrix_two){
  pheno = MotrpacHumanPreSuspensionData::pheno[["data"]] %>%
    dplyr::filter(visitcode == "ADU_BAS") %>%
    dplyr::filter(randomGroupCode == !!randomGroupCode)

  final_output = list()
  shared_cols = intersect(colnames(scion_matrix_one), colnames(scion_matrix_two))
  if(length(shared_cols) == 0 || any(is.na(shared_cols))) stop("desired targets or regulators have no shared participants")

  scion_matrix_one = scion_matrix_one %>%
    dplyr::select(dplyr::all_of(shared_cols)) %>%
    dplyr::select(dplyr::starts_with(as.character(pheno$pid))) #subset to the randomGroupCode

  scion_matrix_two = scion_matrix_two %>%
    dplyr::select(dplyr::all_of(shared_cols)) %>%
    dplyr::select(dplyr::starts_with(as.character(pheno$pid)))  #subset to the randomGroupCode

  if(ncol(scion_matrix_one) == 0 || ncol(scion_matrix_two) == 0)
    stop("desired targets or regulators have no shared participants in this randomGroupCode")

  message(paste("The new matrixes have", ncol(scion_matrix_one),  "shared samples"))
  final_output[["scion_matrix_one"]] = scion_matrix_one
  final_output[["scion_matrix_two"]] = scion_matrix_two
  return(final_output)
}

#' @title Clustering for the scion function
#'
#' @description Scales data and runs clustering to find groups to form networks
#'    from for the scion function
#'
#'
#' @param DA_Input_Cmeans list of differential results that fits the cmeans function
#' @param num_clusters numeric; number of desired clusters.
#'
#' @return A list with relevant clustering result information
#' @export scion_run_cmeans
#'
#' @importFrom Mfuzz mestimate mfuzz
#' @importFrom Biobase ExpressionSet
scion_run_cmeans = function(DA_Input_Cmeans,
                            num_clusters = 13){
  on.exit(gc())

  DA_list <- .prepare_DA_results(DA_list = DA_Input_Cmeans,
                                 convert_features = FALSE,
                                 .contrast_type = "exercise_with_controls")

  nm <- names(DA_list)
  #this part just formats and removes during
  DA_list <- lapply(names(DA_list), function(name_i) {
    zi <- DA_list[[name_i]]

    zi <- zi[, !grepl("during", colnames(zi))]
    zi <- zi[, !grepl("post_10_min", colnames(zi))]

    ome_i <- sub("^[^.]+\\.", "", name_i)
    rownames(zi) <- paste(ome_i, rownames(zi))

    return(zi)
  })

  names(DA_list) <- nm
  #so here is one big difference from the "larger" cmeans function
  #we don't split by tissue and just rbind everything as desired
  zmat_stacked = do.call(rbind, DA_list)

  #same eSet concept, just 1 column of zeros because we're just using 1 modality
  zero_mat <- matrix(data = 0, nrow = nrow(zmat_stacked), ncol = 1L)
  colnames(zero_mat) <- "pre_exercise"
  tmp <- cbind(zero_mat, zmat_stacked)
  # Scale features
  sd <- apply(tmp, 1, sd)
  zmat_stacked <- sweep(zmat_stacked, 1, sd, FUN = "/")
  zmat_stacked <- zmat_stacked[order(sd, decreasing = TRUE), ]
  eSet = Biobase::ExpressionSet(assayData = zmat_stacked)

  #This is still the same stuff as the other c-means stuff, we just call
  #the main mfuzz function from the Mfuzz package.
  m_i <- Mfuzz::mestimate(eSet)

  set.seed(0)
  FCM_i <- Mfuzz::mfuzz(
    eset = eSet,
    centers = num_clusters,
    m = m_i
  )
  return(FCM_i)
}
