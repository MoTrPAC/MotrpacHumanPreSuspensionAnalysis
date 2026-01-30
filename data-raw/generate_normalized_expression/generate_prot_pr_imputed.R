#' Generate imputed protein ratio proteomics QC-normalized datasets
#'
#' This function performs multiple imputation on QC-normalized protein ratio
#' proteomics (Prot-PR) data across all supported tissues. The resulting imputed
#' datasets were used **exclusively for SCION analyses** in the accompanying
#' manuscript and are not used for primary differential analysis workflows.
#'
#' Imputation is carried out using the \code{mice} framework, which supports
#' multiple imputed datasets under a missing-at-random assumption. Parallel
#' execution is supported and strongly recommended due to the computational
#' cost of repeated imputations.
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
#'   \item Loads QC-normalized Prot-PR data
#'   \item Sanitizes column names to satisfy \code{mice} input requirements
#'   \item Performs multiple imputation using \code{run_mice}
#'   \item Restores original sample identifiers
#'   \item Writes imputed matrices to disk with standardized naming
#' }
#'
#' Imputed datasets are written to a dedicated location and clearly labeled to
#' prevent accidental use in standard differential analysis pipelines.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' Imputed protein ratio proteomics matrices are written to disk.
#'
#' @note
#' Imputed Prot-PR datasets are **not** used for the main differential analysis
#' results reported in the manuscript and should only be used for SCION-based
#' analyses.
#'
#' @keywords internal proteomics imputation
#' @author christopher jin

generate_prot_pr_imputed = function(repo_local_dir,
                                    num_cores = 1,
                                    num_imps = 15,
                                    gsutil="gsutil"){
  check_package_installation(pkg = "mice")
  desired_ome = 'prot-pr'; tissue_types = c('muscle', 'adipose')
  output_path = file.path(repo_local_dir, "data", "tmp", "freeze", "proteomics", "qc-norm/")

  for (tissue in tissue_types){
    message(paste("Generating imputed matrixes for prot-pr", tissue))
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
    prot_output_imputed <- run_mice(ome_data,
                                    num_imps=num_imps,
                                    num_cores=num_cores) #by default max allowed missingness is 0.4
    #undo the column name conversion
    colnames(prot_output_imputed) <- gsub("X","",colnames(prot_output_imputed))
    #add ids
    prot_output_imputed <- prot_output_imputed %>%
      dplyr::mutate(feature_id=rownames(prot_output_imputed)) %>%
      dplyr::select(feature_id,everything())
    #save imputed matrix
    write_with_path_name(prot_output_imputed,
                         local_path = output_path,
                         ome = desired_ome,
                         tissue = tissue,
                         data_category = 'imputed',
                         data_details = 'log2-mn')
  }
}
