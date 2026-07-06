#' Prepare and run transcriptomics differential analysis inputs
#'
#' This internal helper function prepares RNA-seq transcriptomics data for
#' differential analysis using linear mixed models. Raw gene-level count data
#' are loaded from cloud storage, filtered to match QC-passed features and
#' samples, and processed using standard \code{edgeR} normalization prior to
#' voom-based model fitting.
#'
#' The function enforces consistency between raw counts and QC-filtered
#' expression data by restricting analysis to genes and participants that
#' pass expression-based filtering thresholds defined upstream.
#'
#' Although both acute and training analyses are supported internally,
#' **only acute exercise transcriptomics results are publicly released** due
#' to sample size and statistical power considerations in training cohorts.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory used for temporary
#' storage and downstream output.
#'
#' @param model_type
#' Character scalar indicating whether the analysis corresponds to
#' \code{"acute"} or \code{"training"} exercise.
#'
#' @param tissue
#' Character scalar specifying the tissue to analyze.
#'
#' @param parallel
#' Logical indicating whether downstream model fitting should be parallelized.
#' Default is \code{FALSE}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Locate and download raw RNA-seq gene count files
#'   \item Load QC-filtered and normalized transcriptomics metadata
#'   \item Restrict genes and samples to those passing expression filters
#'   \item Create an \code{edgeR::DGEList} and compute normalization factors
#'   \item Process covariates and construct mixed model formulas
#'   \item Fit linear mixed models using voom-based precision weighting
#' }
#'
#' Only samples from visit \code{"ADU_BAS"} are used for acute analyses.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' Differential analysis results are written to disk by \code{.run_models}.
#'
#' @note
#' Raw count processing relies on \code{edgeR}. The package is checked at runtime
#' but is not declared as a formal dependency, as transcriptomics analyses are
#' expected to be run by a limited number of analysts.
#'
#' @keywords internal transcriptomics
#' @author christopher jin

.generate_transcriptomics_inputs = function(repo_local_dir,
                                            model_type,
                                            tissue,
                                            parallel = F){
  # check_package_installation("edgeR")
  desired_ome = 'transcript-rna-seq'
  local_path = paste0(repo_local_dir, "data/tmp/")
  counts_data_path = .find_path_name(desired_ome = desired_ome, tissue = tissue, data_type = "rsem-genes-count")
  if(length(counts_data_path) == 0) stop(paste("The desired", tissue, "is not available for transcriptomics"))

  raw_counts_input = MotrpacBicQC::dl_read_gcp(counts_data_path, sep = '\t', tmpdir = local_path)
  parsed_qc_norm = load_qc(selected_omes = desired_ome,
                           selected_tissues = tissue,
                           load_acute_only = FALSE)
  #-> so here we want to make sure the genes we chose and the participants we chose are the same as the ones that pass through expression logcpm cutoffs
  metadata = parsed_qc_norm[[tissue]][[desired_ome]][['sample_metadata']]
  if (model_type == "acute") metadata = metadata %>% dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout
  rownames(metadata) = metadata$vialLabel

  raw_counts_input = raw_counts_input %>%
    dplyr::filter(gene_id %in% rownames(parsed_qc_norm[[tissue]][[desired_ome]][['qc_norm']]))  %>% #filter to select genes
    tibble::column_to_rownames("gene_id") %>% #set to rownames
    dplyr::select(as.character(metadata$vialLabel)) #reorganize and filter to only desired participants

  pre_dge <- edgeR::DGEList(counts = raw_counts_input) #standard dge processing
  pre_dge <- edgeR::calcNormFactors(pre_dge) #standard dge processing
  process_metadata = process_covariates(meta = metadata,
                                        selected_ome = desired_ome,
                                        tissue_input = tissue)

  .run_models(repo_local_dir = repo_local_dir,
              model_type = model_type,
              expression_object = pre_dge,
              process_metadata = process_metadata,
              tissue = tissue,
              ome = desired_ome,
              voom = TRUE,
              parallel = parallel)

}
