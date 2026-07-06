#' Prepare and run acute ATAC-seq differential models for a given tissue
#'
#' This internal helper function prepares ATAC-seq count data and associated
#' covariates for acute exercise analyses and dispatches the appropriate
#' modeling workflow. The function is specific to epigenomic ATAC-seq data
#' generated within the MoTrPAC framework and is intentionally constrained
#' to **acute-only** analyses due to limited availability of training samples.
#'
#' For acute analyses, raw ATAC-seq peak counts are downloaded from cloud
#' storage, filtered to autosomal regions, and aligned to the subset of
#' samples and peaks that pass quality-control and normalization filters.
#' Counts are then converted into an \code{edgeR::DGEList} object with
#' normalization factors computed prior to model fitting.
#'
#' Modeling is performed via \code{.run_models} using a voom-based workflow.
#' No training analysis is performed, even if \code{model_type} is set
#' otherwise, as all available training samples are controls.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory used for temporary
#' storage and downstream processing.
#'
#' @param model_type
#' Character scalar indicating whether the analysis is intended for acute or
#' training comparisons. Only \code{"acute"} is supported for ATAC-seq data.
#'
#' @param tissue
#' Character scalar specifying the tissue to analyze.
#'
#' @param parallel
#' Logical indicating whether model fitting should be parallelized.
#' Default is \code{FALSE}.
#'
#' @details
#' The function performs the following steps for acute analyses:
#' \enumerate{
#'   \item Locate and download raw ATAC-seq peak counts from cloud storage
#'   \item Construct peak identifiers in \code{chr:start-end} format
#'   \item Remove sex chromosome peaks (chrX, chrY)
#'   \item Restrict peaks and samples to those retained after QC normalization
#'   \item Create an \code{edgeR::DGEList} and compute normalization factors
#'   \item Process covariates using \code{process_covariates}
#'   \item Fit voom-based linear models via \code{.run_models}
#' }
#'
#' Only samples from visit \code{"ADU_BAS"} are used, corresponding to the
#' baseline acute exercise bout.
#'
#' @assumptions
#' \itemize{
#'   \item ATAC-seq data are stored under the ome label \code{"epigen-atac-seq"}
#'   \item Peak coordinates are provided as \code{chrom}, \code{start}, and \code{end}
#'   \item Sample identifiers in count matrices correspond to \code{vialLabel}
#'   \item QC-filtered features in \code{qc_norm} define the universe of analyzable peaks
#' }
#'
#' @return
#' This function is called for its side effects and does not return an object.
#' Model outputs are written to disk by \code{.run_models}.
#'
#'
#' @keywords internal
#' @author christopher jin

.generate_atac_inputs = function(repo_local_dir,
                                 model_type,
                                 tissue,
                                 parallel = FALSE){
  check_package_installation("edgeR")
  message("Note: Epigen ATAC Seq only has 5 training samples, all of which are control. No training analysis will be done.")
  desired_ome = 'epigen-atac-seq'
  local_path =  file.path(repo_local_dir, "data/tmp/")
  data_path = .find_path_name(desired_ome = desired_ome, tissue = tissue, data_type = "epigen-atac-seq_counts")
  if (length(data_path) > 0 & model_type == "acute") {
    raw_counts_input = MotrpacBicQC::dl_read_gcp(data_path, sep = '\t', tmpdir = local_path) %>%
      dplyr::mutate(rownames = stringr::str_c(chrom,":",start,"-",end,sep = "")) %>%
      tibble::column_to_rownames("rownames") %>%
      dplyr::filter(!(chrom %in% c("chrX","chrY"))) %>%
      dplyr::select(-chrom,-start,-end) %>%
      as.data.frame()

    parsed_qc_norm = load_qc(selected_omes = desired_ome,
                             selected_tissues = tissue,
                             epigen = TRUE,
                             load_acute_only = FALSE)
    #-> so here we want to make sure the genes we chose and the participants we chose are the same as the ones that pass through expression logcpm cutoffs
    metadata = parsed_qc_norm[[tissue]][[desired_ome]][['sample_metadata']] %>%
      dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout

    rownames(metadata) = metadata$vialLabel
    raw_counts_input = raw_counts_input %>%
      dplyr::filter(rownames(.) %in% rownames(parsed_qc_norm[[tissue]][[desired_ome]][['qc_norm']]))  %>% #filter to select genes
      dplyr::select(as.character(metadata$vialLabel)) #reorganize and filter to only desired participants

    pre_dge <- edgeR::DGEList(counts = raw_counts_input)
    pre_dge <- edgeR::calcNormFactors(pre_dge)

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
}
