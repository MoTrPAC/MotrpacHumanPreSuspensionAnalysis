#' @title Covariates used for statistical analysis
#'
#' @description Table of covariates used in the differential abundance analysis.
#' \code{visitcode} is only used for training analysis. This is also used for
#' any batch correction for regressing out technical covariates in the structure.
#'
#' @usage COVARIATES_FILE
#'
#' @format A \code{data.frame} object with the columns:
#'
#' \describe{
#'   \item{ome}{character; the ome that the covariates are for}
#'   \item{tech_or_design}{factor; either technical or Design, to designate if
#'   individual covariates should be regressed out in the norm-qc tables}
#'   \item{data_type}{factor; which data type a covariate is (e.g. factor)}
#'   \item{covariate}{character; the name of the covariate}
#'   \item{tissue}{factor; one of \code{tissue_available_list}.}
#' }
#'
#' @keywords datasets
"COVARIATES_FILE"
