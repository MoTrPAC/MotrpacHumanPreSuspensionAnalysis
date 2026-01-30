#' Prepare and run metabolomics differential analysis inputs
#'
#' This internal helper function prepares quality-controlled metabolomics data
#' for differential analysis and dispatches model fitting using the standard
#' linear mixed model workflow. Unlike count-based omics layers, metabolomics
#' data are modeled directly from normalized intensity values and therefore do
#' not use voom-based precision weighting.
#'
#' Although the function supports both acute and training analyses internally,
#' **only acute exercise metabolomics results are publicly released** due to
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
#' @param assay
#' Character scalar specifying the metabolomics assay (ome) to process.
#'
#' @param parallel
#' Logical indicating whether downstream model fitting should be parallelized.
#' Default is \code{FALSE}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Load QC-filtered and normalized metabolomics data
#'   \item Optionally restrict samples to baseline acute visits (\code{ADU_BAS})
#'   \item Align expression data with sample metadata
#'   \item Process covariates and construct mixed model formulas
#'   \item Fit linear mixed models via \code{.run_models}
#' }
#'
#' Metabolomics features are filtered to remove unnamed analytes prior to modeling.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' Differential analysis results are written to disk by \code{.run_models}.
#'
#' @note
#' Metabolomics data are modeled on normalized intensities and do not require
#' voom-based weighting.
#'
#' @keywords internal metabolomics
#' @author christopher jin

.generate_metabolomics_inputs = function(repo_local_dir,
                                         acute_vs_training,
                                         tissue,
                                         assay,
                                         parallel = F){
  desired_ome = assay
  ome_data = load_qc(repo_local_dir,
                     selected_tissue = tissue,
                     selected_ome = desired_ome,
                     load_acute_only = FALSE,
                     remove_unnamed_metab = TRUE)
  if (length(ome_data) > 0) {
    metadata = ome_data[[tissue]][[desired_ome]][['sample_metadata']]
    if (acute_vs_training == "acute") metadata = metadata %>% dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout
    rownames(metadata) = metadata$vialLabel
    data_matrix = ome_data[[tissue]][[desired_ome]][['qc_norm']] %>%
      dplyr::select(as.character(metadata$vialLabel))

    process_metadata = process_covariates(meta = metadata,
                                          selected_ome = desired_ome,
                                          tissue_input = tissue,
                                          include_technical = TRUE)

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
