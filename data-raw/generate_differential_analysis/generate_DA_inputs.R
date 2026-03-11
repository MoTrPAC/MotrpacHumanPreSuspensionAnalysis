#' Generate differential analysis (DA) input objects for linear mixed models
#'
#' This function prepares and dispatches preprocessing workflows required to
#' generate analysis-ready inputs for differential analysis using linear mixed
#' models across multiple omics layers and tissues. Preprocessing steps are
#' ome-specific (e.g., raw count handling for RNA-seq, mean-centering for
#' proteomics, assay-specific transformations for metabolomics).
#'
#' The primary purpose of this function is to standardize and expose model-ready
#' inputs in a way that facilitates external methodological development, allowing
#' researchers outside the consortium to evaluate alternative statistical models
#' using harmonized inputs.
#'
#' Due to limited sample sizes and group imbalance, **only acute exercise analyses
#' are publicly released**. Training analyses may be generated internally but are
#' not intended for public dissemination.
#'
#' Methylation-based epigenomic assays are explicitly excluded from this interface.
#' These data are analyzed using generalized linear mixed models rather than the
#' linear mixed model framework used for other omics layers.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory upstream of where
#' intermediate and output files will be written.
#'
#' @param selected_omes
#' Character vector of omics layers to process, or \code{"all"} to include all
#' available omes (excluding unsupported assays).
#'
#' @param selected_tissues
#' Character vector of tissues to process, or \code{"all"} to include all
#' available tissues.
#'
#' @param model_type
#' Character scalar indicating whether inputs should be generated for acute or
#' training analyses. Only \code{"acute"} analyses are publicly released.
#'
#' @param epigen
#' Logical indicating whether epigenomic assays should be included. If \code{FALSE},
#' ATAC-seq inputs are excluded.
#'
#' @param parallel
#' Logical indicating whether downstream model fitting should be parallelized.
#'
#' @details
#' The function iterates over selected tissues and omics layers, dispatching
#' ome-specific helper functions that:
#' \itemize{
#'   \item Load and preprocess raw or normalized data
#'   \item Apply QC-based feature and sample filtering
#'   \item Construct model-ready objects
#'   \item Write standardized DA inputs to disk
#' }
#'
#' Methylation (\code{epigen-methylcap-seq}) is always excluded. ATAC-seq
#' (\code{epigen-atac-seq}) is included only when \code{epigen = TRUE}.
#'
#' @return
#' This function is called for its side effects and does not return an object.
#' All generated inputs are written to disk, organized by tissue and omics layer.
#'
#' @note
#' Training analyses are omitted from public release due to limited sample sizes
#' and lack of adequate group representation.
#'
#' @keywords internal
#' @author christopher jin

generate_DA_inputs = function(repo_local_dir = NULL,
                              selected_omes = "all",
                              selected_tissues = "all",
                              model_type = "acute",
                              epigen = TRUE,
                              parallel = FALSE){
  message("NOTE: Methylation analysis is being processed using a Generalized Linear Mixed Model instead of the Linear Mixed Model that most omes are being processed by. Therefore, performing methylation DA Input will not be offered in this function")
  if(all(selected_omes == "all")) {selected_omes =  ome_available_list()}
  if(all(selected_tissues == "all")){selected_tissues = tissue_available_list()}

  if(!all(selected_omes %in% ome_available_list())){message("Invalid ome selection. Try 'ome_available_list()' for a list")}
  if(!all(selected_tissues %in% tissue_available_list())){message("Invalid tissue selection. Try 'tissue_available_list()' for a list")}

  selected_omes = selected_omes[!selected_omes %in% "epigen-methylcap-seq"] #get rid of methylcap--see above message
  if(!epigen) selected_omes = selected_omes[!selected_omes %in% "epigen-atac-seq"] #remove atac if epigenetics are not desired

  for (tissue in selected_tissues){
    if("transcript-rna-seq" %in% selected_omes) .generate_transcriptomics_inputs(repo_local_dir = repo_local_dir, model_type = model_type, tissue = tissue, parallel = parallel)
    if("prot-ol" %in% selected_omes) .generate_prot_ol_inputs(repo_local_dir = repo_local_dir, model_type = model_type, tissue = tissue, parallel = parallel)
    if("prot-ph" %in% selected_omes) .generate_prot_ph_inputs(repo_local_dir = repo_local_dir, model_type = model_type, tissue = tissue, parallel = parallel)
    if("prot-pr" %in% selected_omes) .generate_prot_pr_inputs(repo_local_dir = repo_local_dir, model_type = model_type, tissue = tissue, parallel = parallel)
    for(metab_assay in grep("^metab-", selected_omes, value = TRUE)) .generate_metabolomics_inputs(repo_local_dir = repo_local_dir, model_type = model_type, tissue = tissue, assay = metab_assay, parallel = parallel)
    if("epigen-atac-seq" %in% selected_omes) .generate_atac_inputs(repo_local_dir = repo_local_dir, model_type = model_type, tissue = tissue, parallel = parallel)
  }
}

#' Run linear mixed model differential analysis and write standardized outputs
#'
#' This internal helper function executes the core differential analysis workflow
#' for a single tissue–omics combination. It fits linear mixed models using the
#' \code{dream} framework, optionally applying voom-based precision weighting,
#' and converts model results into a standardized output format for downstream
#' use and dissemination.
#'
#' While the function supports both acute and training analyses internally,
#' **only acute exercise results are publicly released** due to sample size
#' constraints and statistical power considerations.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory used for output.
#'
#' @param model_type
#' Character scalar indicating whether the model corresponds to acute or training
#' analyses.
#'
#' @param expression_object
#' A model-ready expression object (e.g., \code{edgeR::DGEList} or voom-transformed
#' expression matrix).
#'
#' @param process_metadata
#' A list containing processed sample metadata, including the original metadata
#' and the full model formula.
#'
#' @param tissue
#' Character scalar specifying the tissue being analyzed.
#'
#' @param ome
#' Character scalar specifying the omics layer being analyzed.
#'
#' @param voom
#' Logical indicating whether voom-based precision weights should be applied.
#'
#' @param parallel
#' Logical indicating whether model fitting should be parallelized.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Fit linear mixed models using \code{run_dream}
#'   \item Convert model outputs into a standardized DA format
#'   \item Write results to disk using a consistent naming and directory structure
#' }
#'
#' Outputs are written under \code{data/tmp/freeze_DA/} within the repository
#' directory and are labeled by tissue, omics layer, and analysis type.
#'
#' @return
#' This function is called for its side effects and returns no object.
#'
#' @note
#' Only acute analysis outputs are intended for public release. Training outputs
#' may be generated internally but are not distributed due to limited sample sizes.
#'
#' @keywords internal
#' @author christopher jin

.run_models = function(repo_local_dir,
                       model_type,
                       expression_object,
                       process_metadata,
                       tissue,
                       ome,
                       voom = FALSE,
                       parallel = FALSE){

  da_path = file.path(repo_local_dir, "data", "tmp", "freeze_DA/") #path for output
  dir.create(da_path, recursive = TRUE, showWarnings = FALSE)

  fit = run_dream(expression_object = expression_object,
                  model_type = model_type,
                  process_metadata = process_metadata,
                  voom = voom,
                  parallel = parallel)

  if(model_type == "acute") relevant_formula = process_metadata[["full_formula"]]
  if(model_type == "training") relevant_formula = process_metadata[["training_formula"]]
  if(model_type == "sex_differences") relevant_formula = process_metadata[["sex_differences_formula"]]

  write_output = .convert_dream_output(fit,
                                       metadata = process_metadata$original_meta,
                                       tissue = tissue,
                                       formula = relevant_formula,
                                       ome = ome)
  write_with_path_name(write_output,
                       local_path = da_path,
                       ome = ome,
                       tissue = tissue,
                       data_category = 'da',
                       data_details = paste0('dream-', model_type))

}
