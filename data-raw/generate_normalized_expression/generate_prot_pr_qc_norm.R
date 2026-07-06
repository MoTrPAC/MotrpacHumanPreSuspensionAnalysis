#' Generate QC-normalized protein ratio proteomics (Prot-PR) tables
#'
#' This function generates quality-controlled and normalized protein ratio
#' proteomics (Prot-PR) datasets for downstream analysis. The implementation
#' is derived from the untargeted proteomics pipelines developed by
#' Jim Sanford (adipose) and Patrick Hart (muscle), with tissue-specific
#' handling of technical replicates and batch effects.
#'
#' Raw Prot-PR ratio results are loaded from cloud storage, filtered to remove
#' contaminants and duplicated protein identifiers, normalized using a
#' median-based strategy, and corrected for technical batch effects. Additional
#' post-processing is applied for muscle samples to collapse intra-site and
#' inter-site technical replicates.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory where raw data are
#' downloaded and QC-normalized outputs are written.
#'
#' @details
#' For each supported tissue (muscle and adipose), the function performs:
#' \enumerate{
#'   \item Download raw Prot-PR ratio results and vial-level metadata
#'   \item Remove duplicated protein identifiers and known contaminants
#'   \item Construct a GCT object linking expression data with metadata
#'   \item Restrict analysis to acute study participants
#'   \item Apply median normalization
#'   \item Remove predefined sample outliers
#'   \item Define technical and biological covariates via \code{process_covariates}
#'   \item Remove technical batch effects using \code{limma::removeBatchEffect}
#'   \item Collapse technical replicates (muscle only)
#'   \item Filter features based on missingness thresholds
#'   \item Write QC-normalized matrices and metadata to disk
#' }
#'
#' For muscle tissue, intra-site and inter-site technical replicates are averaged
#' according to predefined naming conventions prior to missingness filtering.
#' For adipose tissue, all samples are processed at a single site and no
#' replicate collapsing is required.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' QC-normalized Prot-PR matrices and associated metadata are written to disk.
#'
#' @note
#' This function requires \code{cmapR} at runtime for GCT object handling.
#' The package is checked dynamically and is not declared as a formal dependency.
#'
#' @keywords internal proteomics normalization
#' @author christopher jin

generate_prot_pr_qc_norm = function(repo_local_dir){
  desired_ome = 'prot-pr'; tissue_types = c('muscle', 'adipose')
  local_path = repo_local_dir
  ome_vial_meta = list()
  ome_vial_meta[["muscle"]] = "gs://motrpac-data-hub/human-precovid/results/proteomics-untargeted/t10-muscle/prot-pr/motrpac_human-precovid_t10-muscle_prot-pr_vial-metadata_v1.0-pnbi.txt"
  ome_vial_meta[["adipose"]] = "gs://motrpac-data-hub/human-precovid/results/proteomics-untargeted/t07-adipose/prot-pr/motrpac_human-precovid_t07-adipose_prot-pr_vial-metadata_v1.0.txt"

  pheno_data_parsed = load_pheno(load_acute_only = FALSE)$pheno_data
  metadata_path = paste0(local_path, "freeze/proteomics/metadata/")
  qc_norm_path = paste0(local_path, "freeze/proteomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  for (tissue in tissue_types){
    message(paste("Generating normalized matrixes for prot-pr", tissue))
    file_load = .find_path_name(desired_ome, tissue = tissue, data_type = 'ratio-results', version = "1.0") #find version 1 of the ratio results
    prot = MotrpacBicQC::dl_read_gcp(file_load, tmpdir = local_path)
    prot = prot[!duplicated(prot$protein_id),] #remove duplicates
    prot = prot %>% filter(prot$is_contaminant == "FALSE")

    mat = prot %>% select(matches("^[0-9]")) %>% as.matrix #raw proteomics
    rdesc <- prot %>% dplyr::select(!matches("^[0-9]"))
    rdesc <- rdesc %>% dplyr::select(!starts_with("pool")) #remove pool info

    rdesc$id <- rdesc$protein_id
    rownames(mat) <- rdesc$protein_id

    tmt_metadata = MotrpacBicQC::dl_read_gcp(ome_vial_meta[[tissue]], tmpdir = local_path) %>%
      dplyr::rename(vialLabel= vial_label) %>%
      dplyr::mutate(id = vialLabel) %>%
      dplyr::mutate(vialLabel = gsub("\\.1", "", vialLabel))
    meta_merge = dplyr::left_join(tmt_metadata, pheno_data_parsed, by = 'vialLabel')
    meta_merge = meta_merge[match(colnames(mat),meta_merge$id),]

    if(ncol(mat)!= nrow(meta_merge) | nrow(mat) != nrow(rdesc)){
      stop("Column or row annotations do not match
		     the matrix size.")}

    #this is non-batch effect corrected, non-normalized dataset
    prot_nonnorm <- cmapR::GCT(mat=mat,
                               rdesc=rdesc,
                               cdesc=meta_merge,
                               rid =rownames(mat),
                               cid = colnames(mat))

    # At some point they changed the indicator for study. the column used to be "protocol" apparently. Not sure when this happened. - Apr 2026, Chris.
    prot_nonnorm_nob <- cmapR::subset_gct(prot_nonnorm, cid = which(prot_nonnorm@cdesc$study == "01")) #filter to just SED
    #This is non-batch effect corrected, median-normalized dataset
    prot_mednorm <- .median_mad_norm(prot_nonnorm_nob, mad = FALSE)

    remove = OUTLIERS$vialLabel
    prot_mednorm_no_outliers <- cmapR::subset_gct(prot_mednorm, cid=which(!(prot_mednorm@cdesc$vialLabel %in% remove)))

    #batch correction
    prot_sample_meta = prot_mednorm_no_outliers@cdesc
    process_metadata = process_covariates(meta = prot_sample_meta,
                                          selected_ome = desired_ome,
                                          tissue_input = tissue)
    meta = process_metadata$metadata
    technical_cov = paste(process_metadata[["technical_cov"]]$covariate, collapse = " + ")
    design_cov = paste(process_metadata[["design_cov"]], collapse = " + ")
    message(tissue, ";", desired_ome, ";technical: ", technical_cov, ";design: ", design_cov)

    prot_mednorm_no_outliers@mat = limma::removeBatchEffect(prot_mednorm_no_outliers@mat,
                                                            covariates = model.matrix(as.formula(paste("~ ", technical_cov)), data = meta),
                                                            design = model.matrix(as.formula(paste("~ ", design_cov)), data = meta))
    #this is batch-corrected, median-normalized dataset
    if(tissue == 'muscle'){
      matr <- prot_mednorm_no_outliers@mat %>% as.data.frame()
      #Per patrick: INTRA site sites end in .1 (and the other version of the rep will have the same first few digits)
      #INTER site replicates end in 7 (and the other version of the rep will have the same first few digits)
      intrasite_replicates = prot_mednorm_no_outliers@cdesc$id[grepl(".1$", prot_mednorm_no_outliers@cdesc$id)]
      for(i in intrasite_replicates){
        j <- gsub("\\.1$", "", i)
        matr[i] <- rowMeans(matr[, c(i, j)], na.rm = TRUE)
      }

      #averaging replicates across chemical sites
      intersite_replicates = prot_mednorm_no_outliers@cdesc$id[grepl("7$", prot_mednorm_no_outliers@cdesc$id)]
      for(j in intersite_replicates){
        i <- paste0(stringr::str_sub(j, end = -2), "2", sep = "")
        matr[i] <- rowMeans(matr[, c(i, j)], na.rm = TRUE)
      }
      #removing columns which have been averaged
      matr <- dplyr::select(matr, -contains(".1"))
      matr <- dplyr::select(matr, -all_of(intersite_replicates))
      matr <- as.matrix(matr)

      cdescc = prot_sample_meta %>%
        dplyr::filter(!grepl("\\.1$", id)) %>%
        dplyr::filter(id %in% colnames(matr))
      cdescc <- cdescc[match(colnames(matr),cdescc$id),]

      if(ncol(matr)!=nrow(cdescc)|nrow(mat)!=nrow(rdesc)){
        stop("Column or row annotations do not match
		     the matrix size. (after avg)")}

      prot_ave <- cmapR::GCT(mat=matr, rdesc= data.frame(rdesc), cdesc=data.frame(cdescc), rid =rownames(matr), cid = colnames(matr))
      PR <- prot_ave %>% .remove_na(0.7)
    }else{
      #For Adipose -> all processing was done at one site
      PR = prot_mednorm_no_outliers %>% .remove_na(0.7)
    }
    protein_pr_output = as.data.frame(PR@mat)
    protein_pr_output$feature_id = rownames(protein_pr_output)
    protein_pr_output = protein_pr_output %>% dplyr::select(feature_id, everything()) #re-order lines

    sample_metadata_output = PR@cdesc %>% dplyr::select(c("vialLabel", "tmt_plex", "tmt16_channel"))

    feature_metadata_output = PR@rdesc %>%
      dplyr::rename(feature_id = protein_id) %>%
      dplyr::select(feature_id, everything()) %>%
      dplyr::filter(feature_id %in% rownames(protein_pr_output)) %>%
      .annotate_prot_pr()

    write_with_path_name(sample_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
    write_with_path_name(protein_pr_output, local_path = qc_norm_path, ome = desired_ome, tissue = tissue, data_category = 'qc-norm', data_details = 'log2-mn')
    write_with_path_name(feature_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features', version = "1.4")
  }
}

.median_mad_norm <- function(x,mad=T){
  #Perform median normalization
  #Args:
  #	x: a GCT object
  #	mad: logic indicating to use mad normalization
  #Returns:
  #	A normalized GCT object
  if(mad){
    scale.factor <- mean(apply(x@mat,2,mad,na.rm=T))
    x@mat <- scale(x@mat,center = apply(x@mat,2,median,na.rm=T),
                   scale = apply(x@mat,2,mad,na.rm=T))
    x@mat <- x@mat*scale.factor
  } else{
    x@mat <- scale(x@mat,center = apply(x@mat,2,median,na.rm=T),
                   scale = F)
  }
  return(x)
}

.remove_na <- function(x,pct){
  #Remove proteins that have more than the set percentage of missing values
  #Args:
  #	x: a GCT object
  #	pct: the percent cutoff
  #Returns:
  #	A GCT object without proteins that have any missing values

  x <- cmapR::subset_gct(x,rid = which(rowSums(is.na(x@mat)) <= pct*ncol(x@mat)))
  return(x)
}

.annotate_prot_pr = function(feature_metadata_output){

  #feature_id is the protein_id (UniProt accession); use it directly from the file as the
  #uniprot key. Strip isoform suffixes (e.g. P12345-2) for the BioMart match.
  prot_pr = feature_metadata_output %>%
    dplyr::mutate(platform = "prot-pr",
                  uniprot = feature_id,
                  uniprot_lookup = stringr::str_remove(uniprot, "-.*"))

  ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprotswissprot")

  prot_lookup_df = biomaRt::getBM(attributes = attributes,
                                   filters = "uniprotswissprot",
                                   values = prot_pr$uniprot_lookup,
                                   mart = ensembl) %>%
    dplyr::distinct() %>%
    dplyr::full_join(prot_pr, by = c("uniprotswissprot" = "uniprot_lookup")) %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(external_gene_name == "", NA, external_gene_name)) %>%
    dplyr::rename(ensembl_gene = ensembl_gene_id,
                  entrez_gene = entrezgene_id) %>%
    #and rearrange to make it ready for the human feature to gene file.
    dplyr::select(assay = platform, feature_id, entrez_gene, gene_symbol, ensembl_gene, uniprot) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::slice_min(entrez_gene, n = 1, with_ties = FALSE)

  #these are a few incorrect gene names that are manually corrected. Not sure how Dan singled out these specific features.
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "SHAN3"] <- "SHANK3"
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "HECD4"] <- "HECTD4"
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "WASH6"] <- "WASH6P"

  return(prot_lookup_df)
}
