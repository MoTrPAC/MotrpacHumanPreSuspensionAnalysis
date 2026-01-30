#' Placeholder for methylation (MethylCap-seq) QC and normalization
#'
#' MethylCap-seq data are **not processed within R** and therefore do not have
#' an associated QC or normalization function in this package. All preprocessing,
#' normalization, and quality-control steps for methylation data are performed
#' using an external, assay-specific pipeline.
#'
#' Users interested in methylation data generation and preprocessing should
#' refer to the Methods section of the accompanying manuscript and consortium
#' documentation for a detailed description of the MethylCap-seq workflow.
#'
#' This function exists solely as a placeholder to maintain a consistent API
#' across omics layers.
#'
#' @return
#' None. This function performs no computation.

generate_methylcap_qc_norm = function(){
}
