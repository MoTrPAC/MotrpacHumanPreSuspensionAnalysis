#' Prepare and run phosphoproteomics differential analysis inputs
#'
#' This internal helper function prepares phosphoproteomics (Prot-PH) data for
#' differential analysis using linear mixed models. Due to non-negligible
#' missingness in phosphoproteomics measurements, features are first filtered
#' using a paired-sample criterion to ensure adequate representation across
#' groups and timepoints.
#'
#' Model fitting is performed using the \code{dream} framework. Because missing
#' data are present, users running DREAM-based analyses must use **R version
#' 4.3 or higher**, which includes the necessary fixes for handling missing values
#' in mixed-model workflows.
#'
#' Although both acute and training analyses are supported internally, **only
#' acute exercise phosphoproteomics results are publicly released** due to
#' limited sample sizes in training cohorts. Certain tissue–ome combinations
#' (e.g., adipose Prot-PH) do not have training samples and are automatically
#' skipped.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory used for data access
#' and output.
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
#'   \item Load QC-filtered phosphoproteomics data
#'   \item Apply paired-sample filtering via \code{filter_paired_n}
#'   \item Optionally restrict samples to baseline acute visits (\code{ADU_BAS})
#'   \item Align expression data with sample metadata
#'   \item Process covariates with technical covariates excluded
#'   \item Fit linear mixed models via \code{.run_models}
#' }
#'
#' Technical covariates are excluded because assay-level normalization already
#' accounts for major sources of technical variation.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' Differential analysis results are written to disk by \code{.run_models}.
#'
#' @note
#' Updated June 2, 2025, to account for reduced sample sizes and missingness
#' characteristics in phosphoproteomics data.
#'
#' @keywords internal phosphoproteomics
#' @author christopher jin

.generate_prot_ph_inputs = function(repo_local_dir,
                                    model_type,
                                    tissue,
                                    parallel = F){
  if(model_type == "training" & tissue == "adipose"){
    message("Note: Adipose Prot-ph/pr has no training samples. No training analysis will be done.")
  }else{
    desired_ome = 'prot-ph'
    ome_data = load_qc(selected_tissues = tissue,
                       selected_omes = desired_ome,
                       load_acute_only = FALSE)
    ome_data = filter_paired_n(qc_data = ome_data,
                               tissue = tissue,
                               ome = desired_ome)

    if (length(ome_data) > 0) {
      metadata = ome_data[[tissue]][[desired_ome]][['sample_metadata']]
      if (model_type == "acute") metadata = metadata %>% dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout

      rownames(metadata) = metadata$vialLabel
      data_matrix = ome_data[[tissue]][[desired_ome]][['qc_norm']] %>%
        dplyr::select(as.character(metadata$vialLabel))

      process_metadata = process_covariates(meta = metadata,
                                            selected_ome = desired_ome,
                                            tissue_input = tissue,
                                            include_technical = F)

      .run_models(repo_local_dir = repo_local_dir,
                  model_type = model_type,
                  expression_object = data_matrix,
                  process_metadata = process_metadata,
                  tissue = tissue,
                  ome = desired_ome,
                  voom = FALSE,
                  parallel = parallel)
    }
  }
}
