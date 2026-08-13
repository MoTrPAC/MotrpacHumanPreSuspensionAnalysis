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
#' which is how the \code{*_DA} objects are labelled.
#'
#' The research platforms are not one object each. They are stacked into a single
#' \code{{TISSUE}_METAB_SUM_STATS} per tissue, which is how
#' \code{{TISSUE}_METAB_DA} is keyed, so the two tiers nest the same way.
#' \code{BLOOD_METAB_T_CLINICAL_SUM_STATS} is clinical chemistry: it carries the
#' same labelling but stays its own object, on this tier and the
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
#' Each object is a \code{data.frame} with one row per feature per randomization
#' group and timepoint, carrying the columns below in this order. They are named
#' and ordered to agree with the \code{*_DA} objects: what the rows are, then what
#' identifies a row, then the statistics.
#'
#' \describe{
#'   \item{tissue}{adipose, blood or muscle.}
#'   \item{assay}{The ome. Every metabolomics platform is stacked under
#'     \code{"metab"} and named in \code{platform}, exactly as the \code{*_DA}
#'     objects do it; every other ome names itself here.}
#'   \item{platform}{The metabolomics platform, for example
#'     \code{"metab-u-rppos"} or \code{"metab-t-clinical"}. Present on
#'     metabolomics objects only, which is also how the \code{*_DA} objects carry
#'     it. Objects for any other ome do not have this column.}
#'   \item{randomGroupCode}{Randomization group: ADUControl, ADUEndur, ADUResist.}
#'   \item{Timepoint}{Timepoint within the acute bout.}
#'   \item{feature_id}{Feature identifier, in the same namespace as the
#'     differential analysis for that tissue and ome.}
#'   \item{Count}{Number of samples summarised, after excluding missing values.}
#'   \item{Mean}{Mean across those samples.}
#'   \item{SD}{Standard deviation across those samples, \code{NA} when
#'     \code{Count} is 1.}
#' }
#'
#' Before v2.0.1 the platform was written into \code{assay}, there was no
#' \code{platform} column, and the identifying columns came first with
#' \code{tissue} and \code{assay} last. There was also one object per
#' metabolomics platform rather than one per tissue.
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
