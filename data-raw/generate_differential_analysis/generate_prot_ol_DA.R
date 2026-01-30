#' Prepare and run Olink proteomics differential analysis inputs
#'
#' This internal helper function prepares Olink proteomics data for differential
#' analysis and dispatches model fitting using the standard linear mixed model
#' workflow. Data are modeled directly from QC-filtered and normalized protein
#' expression values and do not require voom-based precision weighting.
#'
#' Although the function supports both acute and training analyses internally,
#' **only acute exercise proteomics results are publicly released** due to
#' limited sample sizes and statistical power considerations in training cohorts.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory used for data access
#' and output.
#'
#' @param acute_vs_training
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
#'   \item Load QC-filtered and normalized Olink proteomics data
#'   \item Optionally restrict samples to baseline acute visits (\code{ADU_BAS})
#'   \item Align expression data with sample metadata
#'   \item Process covariates with technical covariates excluded
#'   \item Fit linear mixed models via \code{.run_models}
#' }
#'
#' Technical covariates are excluded for Olink proteomics, as assay-specific
#' normalization procedures already account for major technical variation.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' Differential analysis results are written to disk by \code{.run_models}.
#'
#' @note
#' Olink proteomics data are modeled on normalized expression values and do not
#' use voom-based weighting.
#'
#' @keywords internal proteomics
#' @author christopher jin

.generate_prot_ol_inputs = function(repo_local_dir,
                                    acute_vs_training,
                                    tissue,
                                    parallel = F){
  desired_ome = 'prot-ol'
  da_path <- paste0(repo_local_dir, "data/tmp/freeze_DA/") #path for output
  dir.create(da_path, recursive = TRUE, showWarnings = FALSE)

  prot_ol_data = load_qc(repo_local_dir,
                         selected_tissue = tissue,
                         selected_ome = desired_ome,
                         load_acute_only = FALSE)
  if (length(prot_ol_data) > 0) {
    metadata = prot_ol_data[[tissue]][[desired_ome]][['sample_metadata']]
    if (acute_vs_training == "acute") metadata = metadata %>% dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout
    rownames(metadata) = metadata$vialLabel

    data_matrix = prot_ol_data[[tissue]][[desired_ome]][['qc_norm']] %>%
      dplyr::select(as.character(metadata$vialLabel))

    process_metadata = process_covariates(meta = metadata,
                                          selected_ome = desired_ome,
                                          tissue_input = tissue,
                                          include_technical = F)

    .run_models(repo_local_dir = repo_local_dir,
                acute_vs_training = acute_vs_training,
                expression_object = data_matrix,
                process_metadata = process_metadata,
                tissue = tissue,
                ome = desired_ome,
                voom = FALSE,
                parallel = parallel)
  }
}
