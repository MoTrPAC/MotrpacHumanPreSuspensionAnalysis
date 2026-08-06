#' @title Differential Analysis Results
#'
#' @description Differential analysis results. Not to be used directly! Instead,
#'   use \code{\link{load_differential_analysis}}.
#'
#' @usage
#' ADIPOSE_METAB_DA
#' ADIPOSE_PROT_PH_DA
#' ADIPOSE_PROT_PR_DA
#' ADIPOSE_TRNSCRPT_DA
#'
#' BLOOD_METAB_DA
#' BLOOD_PROT_OL_DA
#' BLOOD_METAB_T_CLINICAL_DA
#' BLOOD_PROT_CLINICAL_DA
#' BLOOD_TRNSCRPT_DA
#'
#' MUSCLE_METAB_DA
#' MUSCLE_PROT_PH_DA
#' MUSCLE_PROT_PR_DA
#' MUSCLE_TRNSCRPT_DA
#'
#' @format These objects are not to be used directly. See
#'   \code{\link{load_differential_analysis}} for format details.
#'
#' @keywords datasets
#'
#' @name DA_RESULTS
"ADIPOSE_METAB_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_PROT_PH_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_PROT_PR_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_TRNSCRPT_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_PROT_OL_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_TRNSCRPT_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_PROT_PH_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_PROT_PR_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_TRNSCRPT_DA"



#' Differential Alternative Splicing Results
#'
#' A data object containing results from a differential alternative splicing
#' analysis. The analytical workflow used to generate these results—including
#' read alignment, isoform quantification, splicing event detection, and
#' statistical testing—is described in detail in: "Characterization of exercise-modulated
#' alternative splicing landscape in human skeletal muscle, adipose tissue, and blood in the MoTrPAC Study"

#' To accommodate file size and distribution constraints, this object includes
#' only statistically significant splicing isoforms, defined by a false discovery
#' rate (FDR) less than 0.05. As such, the object represents a filtered subset of
#' all tested isoforms and is intended for downstream analysis, visualization,
#' and integration with other molecular datasets rather than for reproducing the
#' full splicing analysis pipeline.
#'
#' @format A data frame containing differential splicing results for significant
#' isoforms only (FDR < 0.05).
#'
#' @source Generated as described in the splicing companion paper
#'
#' @docType data
#' @keywords datasets
"SPLICING_DA"


#' Clinical Chemistry Differential analysis
#'
#' A data object containing results from a differential expression analysis of the clinical
#' analytes as described in data-raw/differential_analysis_results/
#'
#' @format A \code{data.frame} with 297 rows and 20 columns.
#'
#' @docType data
"CLIN_CHEMISTRY_DA"




#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_T_CLINICAL_DA"

#' @rdname DA_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_PROT_CLINICAL_DA"
