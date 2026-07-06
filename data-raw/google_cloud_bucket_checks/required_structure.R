# Declarative manifest of every file that must be present in the staging bucket.
# Sourced by 05_validate_structure.R.
#
# Columns:
#   ome            — ome string as used in BIC file naming
#   tissue         — high-level tissue category (adipose, blood, muscle)
#   data_category  — one of: qc-norm, metadata, da
#   data_details   — subcategory string in the file name
#   required       — TRUE = pipeline fails if absent; FALSE = warning only
#   gcs_subdir     — top-level GCS subdirectory under PRODUCTION_BUCKET
#   affected_pkg   — "data", "analysis", or "both"
#   tissue_code    — BIC tissue code used in file names (e.g. t06-muscle),
#                    carried directly from OME_TISSUE_CODE

.row = function(ome, tissue, data_category, data_details,
                required = TRUE, gcs_subdir, affected_pkg, tissue_code = NA_character_) {
  data.frame(
    ome = ome, tissue = tissue, tissue_code = tissue_code,
    data_category = data_category, data_details = data_details,
    required = required, gcs_subdir = gcs_subdir, affected_pkg = affected_pkg,
    stringsAsFactors = FALSE
  )
}

# Generates the 4 standard file types for a single ome+tissue combination.
.standard_rows = function(ome, tissue, qc_details, gcs_subdir, tissue_code = NA_character_) {
  rbind(
    .row(ome, tissue, "qc-norm", qc_details, gcs_subdir = gcs_subdir, affected_pkg = "data", tissue_code = tissue_code),
    .row(ome, tissue, "metadata", "features", gcs_subdir = gcs_subdir, affected_pkg = "data", tissue_code = tissue_code),
    .row(ome, tissue, "metadata", "samples", gcs_subdir = gcs_subdir, affected_pkg = "data", tissue_code = tissue_code),
    .row(ome, tissue, "metadata", "removed-samples", required = FALSE, gcs_subdir = gcs_subdir, affected_pkg = "analysis", tissue_code = tissue_code)
  )
}

.da_rows = function(ome, tissue, da_details, gcs_subdir, tissue_code = NA_character_) {
  .row(ome, tissue, "da", da_details, gcs_subdir = gcs_subdir, affected_pkg = "analysis", tissue_code = tissue_code)
}

# Start from OME_TISSUE_CODE and annotate each ome+tissue pair with the
# qc_details string, GCS subdirectory, and DA method for that assay.
# metab-meta-reg is excluded (no staging files).
# clinical-chemistry is not in OME_TISSUE_CODE and is appended manually below.
ome_tissue_pairs = MotrpacHumanPreSuspensionAnalysis::OME_TISSUE_CODE %>%
  dplyr::filter(ome != "metab-meta-reg") %>%
  dplyr::mutate(
    qc_details = dplyr::case_when(
      ome == "transcript-rna-seq"    ~ "log-cpm",
      ome %in% c("prot-pr", "prot-ph") ~ "log2-mn",
      ome == "prot-ol"               ~ "log2",
      ome == "epigen-atac-seq"       ~ "log-cpm",
      ome == "epigen-methylcap-seq"  ~ "beta-values",
      grepl("^metab", ome)           ~ "log2"
    ),
    gcs_subdir = dplyr::case_when(
      ome == "transcript-rna-seq"  ~ "transcriptomics",
      grepl("^prot", ome)          ~ "proteomics",
      grepl("^epigen", ome)        ~ "epigenomics",
      grepl("^metab-u", ome)       ~ "metabolomics-untargeted",
      grepl("^metab-t", ome)       ~ "metabolomics-targeted"
    ),
    da_details = dplyr::case_when(
      ome %in% c("transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol", "epigen-atac-seq") ~ "dream-acute",
      ome == "epigen-methylcap-seq" ~ "malax-glmm-acute"
    )
  )

required_structure = rbind(
  do.call(rbind, lapply(seq_len(nrow(ome_tissue_pairs)), function(i) {
    r = ome_tissue_pairs[i, ]
    rbind(
      .standard_rows(r$ome, r$tissue, r$qc_details, r$gcs_subdir, tissue_code = r$tissue_code),
      if (!is.na(r$da_details)) .da_rows(r$ome, r$tissue, r$da_details, r$gcs_subdir, tissue_code = r$tissue_code) else NULL
    )
  })),

  # clinical-chemistry is not in OME_TISSUE_CODE; tissue_code override is in validate_row()
  .standard_rows("clinical-chemistry", "blood", "log2-transformed", "clinical_chemistry"),
  .da_rows("clinical-chemistry", "blood", "dream-acute", "clinical_chemistry")
)
