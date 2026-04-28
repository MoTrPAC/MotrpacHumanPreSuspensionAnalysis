#' Generate QC-filtered, normalized ATAC-seq matrices for acute analyses
#'
#' This function performs quality control, normalization, and batch correction
#' for human ATAC-seq data and writes standardized outputs to disk for downstream
#' differential analysis. The workflow is specific to acute exercise data and
#' is intentionally restricted to tissues with sufficient sample size.
#'
#' Raw ATAC-seq peak counts are filtered using an aggressive count- and
#' prevalence-based strategy to limit the peak universe to robustly detected
#' regions. Normalization is performed using a voom-based approach with
#' precision weights estimated via the \code{dream} framework, followed by
#' explicit removal of technical batch effects.
#'
#' Training ATAC-seq data are not processed or released due to extremely limited
#' sample sizes (n = 5), which were deemed insufficient for meaningful analysis.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory used for temporary
#' storage and output of normalized matrices and metadata.
#'
#' @param parallel
#' Logical indicating whether voom-based normalization should be parallelized
#' using \code{BiocParallel}. Default is \code{TRUE}.
#'
#' @details
#' The function performs the following steps for each supported tissue:
#' \enumerate{
#'   \item Download raw ATAC-seq peak counts and QA/QC metrics from cloud storage
#'   \item Construct peak identifiers in \code{chr:start-end} format
#'   \item Remove sex chromosome peaks (chrX, chrY)
#'   \item Restrict samples to acute baseline participants
#'   \item Filter peaks based on minimum count and prevalence thresholds
#'   \item Normalize counts using \code{edgeR} and voom precision weighting
#'   \item Remove technical batch effects using \code{limma::removeBatchEffect}
#'   \item Write QC-normalized matrices and metadata to disk
#' }
#'
#' Peak filtering thresholds are defined as:
#' \itemize{
#'   \item Minimum count: \code{2 × median(counts)}
#'   \item Minimum prevalence: \code{50\%} of samples
#' }
#'
#' Covariates used for normalization and batch correction are defined via
#' \code{process_covariates}, with technical covariates removed and biological
#' design variables preserved.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' Normalized matrices and associated metadata are written to disk.
#'
#' @note
#' This function requires \code{variancePartition} and \code{edgeR} at runtime,
#' but these packages are not declared as formal dependencies, as ATAC-seq
#' normalization is expected to be run by a limited number of analysts.
#'
#' @keywords internal atac-seq normalization
#' @author christopher jin


generate_atac_qc_norm = function(repo_local_dir,
                                 parallel = TRUE){
  check_package_installation(pkg = "variancePartition")
  desired_ome = 'epigen-atac-seq'; tissue_types = c('muscle', 'blood')
  local_path = paste0(repo_local_dir, "data/tmp/")
  ome_meta = MotrpacBicQC::dl_read_gcp("gs://motrpac-data-hub/human-precovid/results/epigenomics/qa-qc/motrpac_human-precovid_epigen-atac-seq_qa-qc-metrics_v1.0.csv",
                                       tmpdir = local_path,
                                       sep = ",")
  ome_meta$vialLabel = as.character(ome_meta$vialLabel)

  pheno_data_parsed = load_pheno()$pheno_data
  #so we actually dont load the training here for atac, there's only 5 samples
  #and the group decided that's not usable for anything.
  metadata_path = paste0(local_path, "freeze/epigenomics/metadata/")
  qc_norm_path = paste0(local_path, "freeze/epigenomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  for (tissue in tissue_types){
    message(paste("Generating normalized matrixes for", tissue))
    file_load = .find_path_name(desired_ome, tissue = tissue, data_type = 'count', version = "1.0") #find version 1 of the counts file
    raw_atac_input = MotrpacBicQC::dl_read_gcp(file_load, sep = '\t', tmpdir = local_path)

    raw_atac_input = raw_atac_input %>%
      dplyr::mutate(rownames = stringr::str_c(chrom,":",start,"-",end,sep = "")) %>%
      tibble::column_to_rownames("rownames") %>%
      dplyr::filter(!(chrom %in% c("chrX","chrY"))) %>%
      dplyr::select(-chrom,-start,-end) %>%
      as.data.frame()

    tissue_pheno = pheno_data_parsed[pheno_data_parsed$vialLabel %in% colnames(raw_atac_input), ] #so tissue_pheno has already subsetted only adu-sed participants
    tissue_metadata = ome_meta[ome_meta$vialLabel %in% tissue_pheno$vialLabel, ] #so we use it for the second filtering here
    raw_atac_input = raw_atac_input[, colnames(raw_atac_input) %in% tissue_pheno$vialLabel] #now remove the HA, peds from the counts matrix too

    min_count = 2 * stats::median(as.matrix(raw_atac_input)) #so we go with a more aggresive pruning strategy with ATAC to limit the analyzed peaks
    min_samples = 0.5*dim(raw_atac_input)[2] #number of samples that have to pass the minimum count above

    #---so atac has no outliers, no need to subset those
    atac_filtered = raw_atac_input[rowSums(data.frame(lapply(raw_atac_input, function(x) as.numeric(x >= min_count)), check.names=FALSE)) >= min_samples,]
    meta = merge(tissue_pheno, tissue_metadata, by = "vialLabel") %>%
      tibble::column_to_rownames("vialLabel")
    process_metadata = process_covariates(meta = meta,
                                          selected_ome = desired_ome,
                                          tissue_input = tissue)
    meta = process_metadata$metadata

    dge_list <- edgeR::DGEList(counts = atac_filtered)
    dge_list <- edgeR::calcNormFactors(dge_list)
    formula = process_metadata$full_formula
    message(paste(tissue, "atac", formula))

    if (parallel){
      num_cores = parallel::detectCores() - 2
      param <- BiocParallel::SnowParam(num_cores, "SOCK", progressbar = TRUE)
      suppressWarnings({voom_object <- variancePartition::voomWithDreamWeights(dge_list, formula = stats::as.formula(formula), data = meta, BPPARAM = param)})
    }else{
      voom_object <- variancePartition::voomWithDreamWeights(dge_list, formula = stats::as.formula(formula), data = meta)
    }
    atac_norm <- voom_object$E

    technical_cov = paste(process_metadata[["technical_cov"]]$covariate, collapse = " + ")
    design_cov = paste(process_metadata[["design_cov"]], collapse = " + ")
    message(tissue," technical: ", technical_cov, " design: ", design_cov)
    batch_corrected = limma::removeBatchEffect(atac_norm,
                                               covariates = stats::model.matrix(stats::as.formula(paste("~ ", technical_cov)), data = meta),
                                               design = stats::model.matrix(stats::as.formula(paste("~ ", design_cov)), data = meta))
    batch_corrected = as.data.frame(batch_corrected)
    batch_corrected$feature_id = rownames(batch_corrected)
    batch_corrected = batch_corrected %>%
      dplyr::select(feature_id, everything())
    #for RNA, ATAC, we just want to generate only the list of the features that actually exist so the features in each ome can be easily referenced
    only_features = batch_corrected %>% dplyr::select(feature_id)

    atac_annotated = .annotate_atac_features(only_features)

    write_with_path_name(tissue_metadata, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
    write_with_path_name(batch_corrected, local_path = qc_norm_path, ome = desired_ome, tissue = tissue, data_category = 'qc-norm', data_details = 'log-cpm')
    #in version 1.4 we also attach the gene level information into the feature metadata.
    write_with_path_name(atac_annotated, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features', version = "1.4")

  }
}

.annotate_atac_features = function(feature_metadata){
  atacpeakmeta = feature_metadata %>%
    dplyr::mutate(chrom = gsub(":.*", "", feature_id),
                  start = as.numeric(gsub(".*:|-.*", "", feature_id)),
                  end = as.numeric(gsub(".*-", "", feature_id))) %>%
    data.table::as.data.table()

  atac_peakdf = pre_cawg_get_peak_annotations_hs(atacpeakmeta)
  ensembl = biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 105)

  attributes <- c("ensembl_gene_id", "entrezgene_id", "external_gene_name")
  atac_lookup_df = biomaRt::getBM(attributes = attributes,
                                  filters = "ensembl_gene_id",
                                  values = atac_peakdf$ensembl_gene,
                                  mart = ensembl) %>%
    dplyr::full_join(atac_peakdf, c("ensembl_gene_id" = "ensembl_gene")) %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(external_gene_name == "", NA, external_gene_name)) %>%
    dplyr::select(feature_id,
                  gene_symbol,
                  ensembl_gene = ensembl_gene_id,
                  entrez_gene = entrezgene_id,
                  custom_annotation,
                  relationship_to_gene) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::slice(match(min(entrez_gene), entrez_gene)) %>%
    dplyr::mutate(entrez_gene = as.character(entrez_gene),
                  platform = "epigen-atac-seq")

  return(atac_lookup_df)
}

