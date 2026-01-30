#' @title Contrast Converter
#'
#' @description A \code{data.table} containing information about contrasts.
#'
#' @usage CONTRAST_CONVERTER
#'
#' @format A \code{data.table} with 33 rows and 5 columns, where each row is a
#'   different contrast:
#'
#' \describe{
#'   \item{contrast_order}{integer; the order of the contrast.}
#'   \item{contrast}{factor; full contrast (33 levels).}
#'   \item{contrast_short}{factor; shortened version of the contrasts.}
#'   \item{contrast_type}{factor; one of "exercise_with_controls",
#'   "exercise_no_controls", "Endur_vs_Resist", "baseline", or
#'   "control_only".}
#'   \item{contrast_category}{factor; one of "EE-CON", "RE-CON", "EE-EE",
#'   "RE-RE", "EE-RE", or "CON-CON".}
#' }
#'
#' @keywords datasets
"CONTRAST_CONVERTER"
