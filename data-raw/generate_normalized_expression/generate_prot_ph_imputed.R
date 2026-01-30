#' Generate imputed phosphoproteomics QC-normalized datasets
#'
#' This function performs multiple imputation on QC-normalized phosphoproteomics
#' (Prot-PH) data across all supported tissues. The resulting imputed datasets
#' were used **exclusively for SCION analyses** in the accompanying manuscript
#' and are not used for the primary differential analysis results.
#'
#' Imputation is performed using the \code{mice} framework, which is appropriate
#' for handling missing values under a missing-at-random assumption. Parallel
#' execution is supported and strongly recommended due to the computational
#' cost of multiple imputation.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory where QC-normalized
#' proteomics data are stored and where imputed outputs will be written.
#'
#' @param num_cores
#' Integer specifying the number of cores to use for parallelization during
#' imputation. Default is \code{1}, which disables parallel execution.
#'
#' @param num_imps
#' Integer specifying the number of imputed datasets to generate.
#' Default is \code{15}.
#'
#' @param gsutil
#' Character scalar specifying the path to the \code{gsutil} executable used
#' for data access. Default is \code{"gsutil"}.
#'
#' @details
#' For each tissue, the function:
#' \enumerate{
#'   \item Loads QC-normalized Prot-PH data
#'   \item Temporarily sanitizes column names to meet \code{mice} requirements
#'   \item Performs multiple imputation using \code{run_mice}
#'   \item Restores original sample identifiers
#'   \item Writes the imputed matrix to disk with standardized naming
#' }
#'
#' Imputed values are stored separately from QC-normalized data and are clearly
#' labeled to avoid inadvertent use in primary analyses.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' Imputed phosphoproteomics matrices are written to disk.
#'
#' @note
#' Imputed Prot-PH datasets are **not** used for differential analysis results
#' reported in the main manuscript and should not be substituted for QC-normalized
#' data in standard modeling workflows.
#'
#' @keywords internal proteomics imputation
#' @author christopher jin

generate_prot_ph_imputed = function(repo_local_dir,
                                    num_cores = 1,
                                    num_imps = 15,
                                    gsutil="gsutil"){
  check_package_installation(pkg = "mice")
  desired_ome = 'prot-ph'; tissue_types = c('muscle', 'adipose')
  output_path = file.path(repo_local_dir, "data", "tmp", "freeze", "proteomics", "qc-norm/")

  for (tissue in tissue_types){
    message(paste("Generating imputed matrixes for prot-ph", tissue))
    #get QC norm tables
    data_list = load_qc(repo_local_dir,
                        selected_tissues = tissue,
                        selected_omes = desired_ome,
                        load_acute_only = FALSE,
                        gsutil = gsutil)
    ome_data = data_list[[tissue]][[desired_ome]][["qc_norm"]]
    #mice requires non-numeric column names
    colnames(ome_data) <- make.names(colnames(ome_data))
    #perform imputation
    #parallelization highly recommended (adjust num_cores parameter)
    phospho_output_imputed <- run_mice(ome_data,
                                       num_imps=num_imps,
                                       num_cores=num_cores)
    #undo the column name conversion
    colnames(phospho_output_imputed) <- gsub("X","",colnames(phospho_output_imputed))
    #add ids
    phospho_output_imputed <- phospho_output_imputed %>%
      dplyr::mutate(feature_id=rownames(phospho_output_imputed)) %>%
      dplyr::select(feature_id,everything())
    #save imputed matrix
    write_with_path_name(phospho_output_imputed,
                         local_path = output_path,
                         ome = desired_ome,
                         tissue = tissue,
                         data_category = 'imputed',
                         data_details = 'log2-mn')
  }
}
