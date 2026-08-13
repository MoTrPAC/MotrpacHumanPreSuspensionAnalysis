#' @title Summary Statistics by Tissue, Assay, and Platform
#'
#' @description
#' Group- and timepoint-level summary statistics computed across multiple
#' tissues, molecular assays, and analytical platforms. Each dataset contains
#' means and standard deviations calculated within randomization groups and
#' timepoints. These objects are intended for descriptive and exploratory use
#' and are not differential analysis results. Available for all omes/tissue/platforms
#' except methylcap, which is processed separately.
#'
#' For any datasets with missing values, the values reflect means/sds/n with missing values excluded.
#'
#' @details
#' Datasets are stratified by tissue (e.g., adipose, blood, muscle) and assay
#' (metabolomics, proteomics, transcriptomics, epigenomics) — one object per
#' tissue and assay, named the way the \code{*_DA} objects are named.
#'
#' Every metabolomics object reads \code{assay = "metab"} and names its platform
#' in a \code{platform} column (e.g. \code{"metab-u-rppos"}, \code{"metab-t-tca"}),
#' which is how the \code{*_DA} objects are labelled. Before v2.0.1 this tier
#' named the platform in \code{assay} instead, so the two tiers disagreed about
#' what an ome was called and every join between them had to translate.
#'
#' The research platforms are stacked into a single
#' \code{{TISSUE}_METAB_SUM_STATS} per tissue.
#' \code{BLOOD_METAB_T_CLINICAL_SUM_STATS} is clinical chemistry: it is labelled
#' the same way but stays its own object, on this tier and the
#' differential-analysis tier alike, because it shares five analytes with the
#' research platforms (Cortisol, Glucose, Glycerol, KET, NEFA) that would
#' otherwise sit in one object twice. \code{platform} is what tells the two
#' apart — \code{assay} does not.
#'
#' @usage
#' ## Adipose
#' ADIPOSE_EPIGEN_METHYLCAP_SEQ_SUM_STATS
#' ADIPOSE_METAB_SUM_STATS
#'
#' ADIPOSE_PROT_PH_SUM_STATS
#' ADIPOSE_PROT_PR_SUM_STATS
#' ADIPOSE_TRANSCRIPT_RNA_SEQ_SUM_STATS
#'
#' ## Blood
#' BLOOD_EPIGEN_METHYLCAP_SEQ_SUM_STATS
#' BLOOD_METAB_SUM_STATS
#'
#' BLOOD_PROT_OL_SUM_STATS
#' BLOOD_METAB_T_CLINICAL_SUM_STATS
#' BLOOD_PROT_CLINICAL_SUM_STATS
#' BLOOD_TRANSCRIPT_RNA_SEQ_SUM_STATS
#'
#' ## Muscle
#' MUSCLE_EPIGEN_METHYLCAP_SEQ_SUM_STATS
#' MUSCLE_METAB_SUM_STATS
#'
#' MUSCLE_PROT_PH_SUM_STATS
#' MUSCLE_PROT_PR_SUM_STATS
#' MUSCLE_TRANSCRIPT_RNA_SEQ_SUM_STATS
#' MUSCLE_EPIGEN_ATAC_SEQ_SUM_STATS
#'
#' @format
#' Each object is a \code{data.frame} with one row per feature per
#' randomization group and timepoint, carrying \code{randomGroupCode},
#' \code{feature_id}, \code{Timepoint}, \code{Count}, \code{Mean}, \code{SD},
#' \code{tissue} and \code{assay}. The metabolomics objects — the three
#' \code{{TISSUE}_METAB_SUM_STATS} and \code{BLOOD_METAB_T_CLINICAL_SUM_STATS} —
#' carry one further column, \code{platform}, and order their columns
#' the way the \code{*_DA} objects order the ones they share:
#' \code{tissue}, \code{assay}, \code{platform}, \code{randomGroupCode},
#' \code{Timepoint}, \code{feature_id}, \code{Count}, \code{Mean}, \code{SD}.
#'
#' @keywords datasets
#'
#' @name SUM_STATS_RESULTS
"ADIPOSE_EPIGEN_METHYLCAP_SEQ_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_PROT_PH_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_PROT_PR_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_TRANSCRIPT_RNA_SEQ_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_EPIGEN_METHYLCAP_SEQ_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_PROT_OL_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_TRANSCRIPT_RNA_SEQ_SUM_STATS"


#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_EPIGEN_METHYLCAP_SEQ_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_PROT_PH_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_PROT_PR_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_TRANSCRIPT_RNA_SEQ_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_EPIGEN_ATAC_SEQ_SUM_STATS"


#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_T_CLINICAL_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_PROT_CLINICAL_SUM_STATS"
