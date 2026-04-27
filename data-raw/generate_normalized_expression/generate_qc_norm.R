#' Generate QC-normalized datasets for the human Pre-COVID analysis
#'
#' This wrapper function orchestrates quality control, normalization, and
#' (where applicable) imputation across all supported omics layers for the
#' human Pre-COVID dataset. It sequentially calls ome-specific normalization
#' functions and writes standardized outputs to disk, which together constitute
#' the data objects accessed via \code{\link{load_qc}}.
#'
#' The function is intended to be run once per data freeze or release cycle and
#' provides a reproducible, end-to-end mechanism for generating all analysis-ready
#' inputs used in downstream modeling and visualization workflows.
#'
#' @details
#' Depending on user-specified options, this function may generate:
#' \itemize{
#'   \item QC-normalized transcriptomics data (RNA-seq)
#'   \item QC-normalized proteomics data (Prot-OL, Prot-PH, Prot-PR)
#'   \item Imputed proteomics datasets (Prot-PH, Prot-PR; SCION analyses only)
#'   \item QC-normalized metabolomics data
#'   \item QC-normalized ATAC-seq data (epigenomics)
#' }
#'
#' Imputed datasets are generated **only** for SCION analyses and are not used
#' for primary differential analysis results. Training-phase analyses and
#' normalization steps are intentionally excluded or limited for certain omes
#' due to sample size and power considerations.
#'
#' For details on ome-specific processing steps, refer to the corresponding
#' functions:
#' \describe{
#'   \item{\code{\link{generate_transcriptomics_qc_norm}}}{RNA-seq normalization and filtering}
#'   \item{\code{\link{generate_prot_ol_qc_norm}}}{Olink targeted proteomics normalization}
#'   \item{\code{\link{generate_prot_ph_qc_norm}}}{Phosphoproteomics normalization}
#'   \item{\code{\link{generate_prot_ph_imputed}}}{Phosphoproteomics imputation (SCION only)}
#'   \item{\code{\link{generate_prot_pr_qc_norm}}}{Protein ratio proteomics normalization}
#'   \item{\code{\link{generate_prot_pr_imputed}}}{Protein ratio proteomics imputation (SCION only)}
#'   \item{\code{\link{generate_atac_qc_norm}}}{ATAC-seq normalization}
#' }
#'
#' @param repo_local_dir
#' Character scalar specifying the local directory where all intermediate and
#' final QC-normalized outputs will be written.
#'
#' @param selected_omes
#' Character vector specifying which omics layers to generate. Use
#' \code{\link{ome_available_list}} to view valid options.
#'
#' @param selected_tissues
#' Character vector specifying which tissues to generate. Use
#' \code{\link{tissue_available_list}} to view valid options.
#'
#' @param epigen
#' Logical indicating whether epigenomic datasets (e.g., ATAC-seq) should be
#' generated. Setting this to \code{TRUE} substantially increases compute time.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' Successful completion indicates that all requested QC-normalized datasets
#' have been written to disk and are available via \code{\link{load_qc}}.
#'
#' @author Christopher Jin
#'
#' @importFrom MotrpacHumanPreSuspensionData ome_available_list tissue_available_list

generate_qc_norm = function(repo_local_dir = NULL,
                            selected_omes = "all",
                            selected_tissues = "all",
                            epigen = TRUE,
                            config_file = NULL){
  if(all(selected_omes == "all")) selected_omes = ome_available_list()
  if(all(selected_tissues == "all")) selected_tissues = tissue_available_list()

  if(!all(selected_omes %in% ome_available_list())){
    stop("Invalid ome selection. Try 'ome_available_list()' for a list")
  }

  if(!all(selected_tissues %in% tissue_available_list())){
    stop("Invalid tissue selection. Try 'tissue_available_list()' for a list")
  }

  if(!epigen){
    selected_omes <- selected_omes[!selected_omes %in% "epigen-atac-seq"]
    selected_omes <- selected_omes[!selected_omes %in% "epigen-methylcap-seq"]
  } #remove atac, methyl if epigenetics are not desired

  if("transcript-rna-seq" %in% selected_omes) generate_transcriptomics_qc_norm(repo_local_dir)
  if("prot-ol" %in% selected_omes) generate_prot_ol_qc_norm(repo_local_dir)
  if("prot-ph" %in% selected_omes){
    generate_prot_ph_qc_norm(repo_local_dir)
    generate_prot_ph_imputed(repo_local_dir)
  }
  if("prot-pr" %in% selected_omes){
    generate_prot_pr_qc_norm(repo_local_dir)
    generate_prot_pr_imputed(repo_local_dir)
  }
  if(any(grepl("metab", selected_omes))) {
    #this function requires a 'config' file because it needs to source files from multiple repositories
   generate_metab_qc_norm(config_file = config_file)
  }
  if("epigen-atac-seq" %in% selected_omes) generate_atac_qc_norm(repo_local_dir)
  if("epigen-methyl-seq" %in% selected_omes){
    message("Methylation processing requires analysis pipelines not compatible
            with this generate function. Please reach out to cajin@stanford.edu
            or yongchao.ge@mssm.edu if any questions about methylation come up.")
  }
  return("QC Norm Generation Complete")
}


