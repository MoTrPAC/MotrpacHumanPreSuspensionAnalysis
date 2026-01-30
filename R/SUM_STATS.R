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
#' Datasets are stratified by tissue (e.g., adipose, blood, muscle), assay
#' (metabolomics, proteomics, transcriptomics, epigenomics), and where
#' applicable, analytical platform (e.g., TCA, AMINES, OXYLIPNEG, RPNEG,
#' RPPOS, ATAC-SEQ, RNA-SEQ).
#'
#' @usage
#' ## Adipose
#' ADIPOSE_METAB_T_ACOA_SUM_STATS
#' ADIPOSE_METAB_T_AMINES_SUM_STATS
#' ADIPOSE_METAB_T_TCA_SUM_STATS
#' ADIPOSE_METAB_T_NUC_SUM_STATS
#' ADIPOSE_METAB_T_KA_SUM_STATS
#' ADIPOSE_METAB_T_OXYLIPNEG_SUM_STATS
#' ADIPOSE_METAB_U_HILICPOS_SUM_STATS
#' ADIPOSE_METAB_U_IONPNEG_SUM_STATS
#' ADIPOSE_METAB_U_LRPNEG_SUM_STATS
#' ADIPOSE_METAB_U_LRPPOS_SUM_STATS
#' ADIPOSE_METAB_U_RPNEG_SUM_STATS
#' ADIPOSE_METAB_U_RPPOS_SUM_STATS
#'
#' ADIPOSE_PROT_PH_SUM_STATS
#' ADIPOSE_PROT_PR_SUM_STATS
#' ADIPOSE_TRANSCRIPT_RNA_SEQ_SUM_STATS
#'
#' ## Blood
#' BLOOD_METAB_T_AMINES_SUM_STATS
#' BLOOD_METAB_T_TCA_SUM_STATS
#' BLOOD_METAB_T_CONV_SUM_STATS
#' BLOOD_METAB_T_OXYLIPNEG_SUM_STATS
#' BLOOD_METAB_U_HILICPOS_SUM_STATS
#' BLOOD_METAB_U_IONPNEG_SUM_STATS
#' BLOOD_METAB_U_LRPNEG_SUM_STATS
#' BLOOD_METAB_U_LRPPOS_SUM_STATS
#' BLOOD_METAB_U_RPNEG_SUM_STATS
#' BLOOD_METAB_U_RPPOS_SUM_STATS
#'
#' BLOOD_PROT_OL_SUM_STATS
#' BLOOD_TRANSCRIPT_RNA_SEQ_SUM_STATS
#' BLOOD_EPIGEN_ATAC_SEQ_SUM_STATS
#'
#' ## Muscle
#' MUSCLE_METAB_T_AMINES_SUM_STATS
#' MUSCLE_METAB_T_TCA_SUM_STATS
#' MUSCLE_METAB_T_NUC_SUM_STATS
#' MUSCLE_METAB_T_OXYLIPNEG_SUM_STATS
#' MUSCLE_METAB_U_HILICPOS_SUM_STATS
#' MUSCLE_METAB_U_IONPNEG_SUM_STATS
#' MUSCLE_METAB_U_LRPNEG_SUM_STATS
#' MUSCLE_METAB_U_LRPPOS_SUM_STATS
#' MUSCLE_METAB_U_RPNEG_SUM_STATS
#' MUSCLE_METAB_U_RPPOS_SUM_STATS
#'
#' MUSCLE_PROT_PH_SUM_STATS
#' MUSCLE_PROT_PR_SUM_STATS
#' MUSCLE_TRANSCRIPT_RNA_SEQ_SUM_STATS
#' MUSCLE_EPIGEN_ATAC_SEQ_SUM_STATS
#'
#' @format
#' Each object is a \code{data.frame} with one row per feature per
#' randomization group and timepoint. Columns typically include feature
#' identifiers, tissue, assay, platform, group, timepoint, mean, and standard
#' deviation.
#'
#' @keywords datasets
#'
#' @name SUM_STATS_RESULTS
"ADIPOSE_METAB_T_ACOA_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_T_AMINES_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_T_KA_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_T_NUC_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_T_OXYLIPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_T_TCA_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_U_HILICPOS_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_U_IONPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_U_LRPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_U_LRPPOS_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_U_RPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"ADIPOSE_METAB_U_RPPOS_SUM_STATS"

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
"BLOOD_METAB_T_AMINES_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_T_CONV_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_T_OXYLIPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_T_TCA_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_U_HILICPOS_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_U_IONPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_U_LRPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_U_LRPPOS_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_U_RPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"BLOOD_METAB_U_RPPOS_SUM_STATS"

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
"BLOOD_EPIGEN_ATAC_SEQ_SUM_STATS"



#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_T_AMINES_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_T_NUC_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_T_OXYLIPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_T_TCA_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_U_HILICPOS_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_U_IONPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_U_LRPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_U_LRPPOS_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_U_RPNEG_SUM_STATS"

#' @rdname SUM_STATS_RESULTS
#' @format NULL
#' @usage NULL
"MUSCLE_METAB_U_RPPOS_SUM_STATS"

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

