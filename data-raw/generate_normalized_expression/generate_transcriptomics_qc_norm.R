#' Generate QC-normalized transcriptomics datasets for Pre-COVID analysis
#'
#' This function generates quality-controlled and normalized RNA-seq expression
#' matrices for downstream visualization, clustering, and exploratory analyses
#' in the human Pre-COVID dataset. It is a subroutine of the higher-level
#' \code{generate_qc_norm_data} workflow and produces the transcriptomics
#' components accessed via \code{\link{load_qc}}.
#'
#' Raw gene-level count data are filtered to remove lowly expressed genes,
#' normalized using the trimmed mean of M-values (TMM) method, transformed to
#' log-counts-per-million (log-CPM), and batch-corrected for technical covariates.
#'
#' @details
#' For each supported tissue, the function performs the following steps:
#' \enumerate{
#'   \item Download raw RNA-seq gene count matrices from cloud storage
#'   \item Remove predefined sample outliers
#'   \item Merge sequencing QA/QC metadata with phenotypic metadata
#'   \item Filter genes with low expression using CPM-based thresholds
#'   \item Normalize library sizes using TMM normalization
#'   \item Transform normalized counts to log-CPM
#'   \item Remove technical batch effects while preserving biological design factors
#'   \item Write QC-normalized matrices and associated metadata to disk
#' }
#'
#' Gene filtering criteria are defined as:
#' \itemize{
#'   \item Expression > 0.5 CPM
#'   \item Present in at least 10\% of samples within a tissue
#' }
#'
#' Batch correction is performed using \code{limma::removeBatchEffect} and is
#' intended **only for visualization and exploratory analyses**. Raw count data
#' and uncorrected normalization outputs are used for differential analysis.
#'
#' Outputs are written in a directory structure compatible with direct
#' synchronization (e.g., via \code{rsync}) to the MoTrPAC data hub.
#'
#' @param repo_local_dir
#' Character scalar specifying the local directory where normalized expression
#' matrices and metadata will be written.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' QC-normalized transcriptomics matrices and associated metadata are written
#' to disk.
#'
#' @note
#' Differential analysis workflows should operate on raw counts or voom-weighted
#' inputs generated elsewhere in the pipeline. The batch-corrected log-CPM
#' matrices produced here are not used for statistical testing.
#'
#' @keywords internal transcriptomics normalization
#' @author christopher jin

generate_transcriptomics_qc_norm = function(repo_local_dir){
  desired_ome = 'transcript-rna-seq'
  tissue_types = c('muscle', 'blood', 'adipose')
  local_path = repo_local_dir

  merged_metadata = MotrpacBicQC::dl_read_gcp('gs://motrpac-data-hub/human-precovid/results/transcriptomics/qa-qc/motrpac_human-precovid_transcript-rna-seq_qa-qc-metrics_v1.0.csv',
                                              tmpdir = local_path, sep = ",", check_first = TRUE) %>%
    dplyr::select(-"BID") %>%  #just remove it here to avoid these columns from becoming duplicate column names
    dplyr::select(-"PID") #just remove it here to avoid these columns from becoming duplicate column names

  merged_metadata$Batch <- sub(".*?(\\d{1,2})$", "\\1", merged_metadata$Lib_batch_ID) #choose last 2 characters
  pheno_data_parsed = load_pheno(load_acute_only = FALSE)$pheno_data

  metadata_path = paste0(local_path, "freeze/transcriptomics/metadata/") #path for output
  qc_norm_path = paste0(local_path, "freeze/transcriptomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  for (tissue in tissue_types){
    message(paste("Generating transcriptomics normalized matrixes for", tissue))
    file_load = .find_path_name(desired_ome, tissue = tissue, data_type = 'count', version = "1.0") #find version 1 of the counts file
    raw_counts_input = MotrpacBicQC::dl_read_gcp(file_load, sep = '\t', tmpdir = local_path) %>%
      tibble::column_to_rownames("gene_id")

    outliers_vialLabels = OUTLIERS$vialLabel
    desired_samples = setdiff(colnames(raw_counts_input), outliers_vialLabels) #basically the counts minus the outliers
    tissue_pheno = pheno_data_parsed[pheno_data_parsed$vialLabel %in% desired_samples, ] #so tissue_pheno has already subsetted only adu-sed participants
    tissue_metadata = merged_metadata[merged_metadata$vialLabel %in% tissue_pheno$vialLabel, ] #so we use it for the second filtering here
    raw_counts_input = raw_counts_input[, colnames(raw_counts_input) %in% tissue_pheno$vialLabel] #now remove the HA, peds from the counts matrix too

    meta = merge(tissue_pheno, tissue_metadata, by = "vialLabel")
    process_metadata = process_covariates(meta = meta,
                                          selected_ome = desired_ome,
                                          tissue_input = tissue)
    meta = process_metadata$metadata

    #-----here we filter lowly expressed genes and transform into logcpm
    raw_dge = edgeR::DGEList(counts = raw_counts_input)
    keep = rowSums(cpm(raw_dge) > 0.5) >= round(length(raw_counts_input)*0.1)
    filt_dge = raw_dge[keep, , keep.lib.sizes=FALSE]
    dge = edgeR::calcNormFactors(filt_dge, method="TMM")
    norm_counts = edgeR::cpm(dge,log=TRUE)

    technical_cov = paste(process_metadata[["technical_cov"]]$covariate, collapse = " + ")
    design_cov = paste(process_metadata[["design_cov"]], collapse = " + ")
    message(tissue," technical: ", technical_cov, " design: ", design_cov)
    #---here we perform batch correction:: ONLY for visualization/clustering/etc. Use raw counts for DA
    batch_corrected = limma::removeBatchEffect(norm_counts,
                                               covariates = model.matrix(as.formula(paste("~ ", technical_cov)), data = meta),
                                               design = model.matrix(as.formula(paste("~ ", design_cov)), data = meta))
    batch_corrected = as.data.frame(batch_corrected)
    batch_corrected$feature_id = rownames(batch_corrected)
    batch_corrected = batch_corrected %>%
      dplyr::select(feature_id, everything())

    #for RNA, ATAC, we just want to generate only the list of the features that actually exist so the features in each ome can be easily referenced
    only_features = batch_corrected %>% dplyr::select(feature_id)

    write_with_path_name(tissue_metadata, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
    write_with_path_name(batch_corrected, local_path = qc_norm_path, ome = desired_ome, tissue = tissue, data_category = 'qc-norm', data_details = 'log-cpm')
    write_with_path_name(only_features, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features')
  }
}

