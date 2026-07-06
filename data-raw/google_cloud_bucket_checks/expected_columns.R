# Per-ome expected column schemas for structural validation (Step 5).
# Sourced by 05_validate_structure.R.
#
# Each entry is a named list with:
#   required  — columns that must be present
#   forbidden — columns that must NOT be present (catches stale schema)
#
# DA files all share the same core schema regardless of ome.
# qc-norm and metadata schemas are ome-specific.

DA_COLS_CORE = c(
  "feature_id", "contrast", "full_model",
  "logFC", "CI.L", "CI.R",
  "t", "p_value", "adj_p_value"
)

DA_COLS_PROT = c(DA_COLS_CORE, "AveExpr", "degrees_of_freedom", "logLik")

DA_COLS_METAB = c(DA_COLS_CORE, "AveExpr", "platform")

DA_COLS_METHYL = c(
  "feature_id", "contrast", "full_model",
  "methylation_diff", "p_value", "adj_p_value"
)

QC_NORM_COLS_EXPR = c("feature_id") # + sample ID columns (checked by nrow/ncol)

#since we combine with pheno upon loading, all we need is the vialLabel for some sample mapping.
METADATA_SAMPLE_COLS_CORE = c(
  "vialLabel"
)

METADATA_FEATURE_COLS_TRANSCRIPT = c(
  "assay", "feature_id", "gene_symbol", "ensembl_gene"
)

METADATA_FEATURE_COLS_PROT = c(
  "assay", "feature_id", "gene_symbol", "uniprot"
)

METADATA_FEATURE_COLS_METAB = c(
  "assay", "feature_id", "refmet_name"
)

METADATA_FEATURE_COLS_EPIGEN = c(
  "assay", "feature_id", "gene_symbol", "ensembl_gene", "entrez_gene", "relationship_to_gene"
)

# Map from (data_category, ome) to the set of required columns.
# Keys are "<data_category>__<ome>" or "<data_category>__*" for shared schemas.
EXPECTED_COLUMNS = list(

  # DA schemas
  "da__transcript-rna-seq"  = list(required = DA_COLS_CORE,  forbidden = character(0)),
  "da__prot-pr"             = list(required = DA_COLS_PROT,  forbidden = character(0)),
  "da__prot-ph"             = list(required = DA_COLS_PROT,  forbidden = character(0)),
  "da__prot-ol"             = list(required = DA_COLS_CORE,  forbidden = character(0)),
  "da__metab"               = list(required = DA_COLS_METAB, forbidden = character(0)),
  "da__epigen-methylcap-seq"= list(required = DA_COLS_METHYL,forbidden = character(0)),
  "da__clinical-chemistry"  = list(required = DA_COLS_CORE,  forbidden = character(0)),

  # qc-norm schemas (expression matrix: just need feature_id + at least 1 sample col)
  "qc-norm__transcript-rna-seq"   = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__prot-pr"              = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__prot-ph"              = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__prot-ol"              = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__epigen-atac-seq"      = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__epigen-methylcap-seq" = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__metab"                = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__clinical-chemistry"   = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),

  # metadata schemas
  "metadata__samples"         = list(required = METADATA_SAMPLE_COLS_CORE, forbidden = character(0)),
  # removed-samples files use either "vialLabel" or "sample" as the ID column
  "metadata__removed-samples" = list(required = character(0), forbidden = character(0)),
  "metadata__features__transcript-rna-seq" = list(required = METADATA_FEATURE_COLS_TRANSCRIPT,  forbidden = character(0)),
  "metadata__features__prot-pr"            = list(required = METADATA_FEATURE_COLS_PROT,        forbidden = character(0)),
  "metadata__features__prot-ph"            = list(required = METADATA_FEATURE_COLS_PROT,        forbidden = character(0)),
  "metadata__features__prot-ol"            = list(required = METADATA_FEATURE_COLS_PROT,        forbidden = character(0)),
  "metadata__features__epigen-atac-seq"     = list(required = METADATA_FEATURE_COLS_EPIGEN, forbidden = character(0)),
  "metadata__features__epigen-methylcap-seq" = list(required = METADATA_FEATURE_COLS_EPIGEN, forbidden = character(0)),
  "metadata__features__metab"              = list(required = METADATA_FEATURE_COLS_METAB,       forbidden = character(0))
)

# Helper: look up the expected columns for a given data_category + ome combination.
# Falls back to the generic metab key for any metab platform.
.get_expected_cols = function(data_category, ome, data_details = NULL) {
  key = paste0(data_category, "__", ome)
  if (!is.null(key) && key %in% names(EXPECTED_COLUMNS)) {
    return(EXPECTED_COLUMNS[[key]])
  }
  # metab fallback — try with data_details first (e.g. "metadata__features__metab"),
  # then without (e.g. "qc-norm__metab", "da__metab")
  if (grepl("^metab", ome)) {
    if (data_category == "metadata" && !is.null(data_details)) {
      key_metab_detail = paste0("metadata__", data_details, "__metab")
      if (key_metab_detail %in% names(EXPECTED_COLUMNS)) return(EXPECTED_COLUMNS[[key_metab_detail]])
    }
    key_metab = paste0(data_category, "__metab")
    if (key_metab %in% names(EXPECTED_COLUMNS)) {
      return(EXPECTED_COLUMNS[[key_metab]])
    }
  }
  # metadata samples vs features
  if (data_category == "metadata" && !is.null(data_details)) {
    key_meta = paste0("metadata__", data_details, "__", ome)
    if (key_meta %in% names(EXPECTED_COLUMNS)) return(EXPECTED_COLUMNS[[key_meta]])
    key_meta_shared = paste0("metadata__", data_details)
    if (key_meta_shared %in% names(EXPECTED_COLUMNS)) return(EXPECTED_COLUMNS[[key_meta_shared]])
  }
  return(NULL)
}
