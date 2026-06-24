# Step 5: Validate the staging bucket structure and values against the required manifest.
#
# Single-pass design: each staging file is downloaded exactly once. validate_row()
# checks column schema and, for qc-norm / metadata-features / metadata-samples rows,
# caches the full data frame in staging_data_cache. After all rows are validated the
# value comparison loop reads from that cache — no second download.
#
# Structure checks:
#   Checks for required ome x tissue x data_category combinations and column schemas.
#   Selects the highest-versioned file when multiple versions are present.
#   Writes logs/structure_validation_{timestamp}.tsv.
#   Pipeline fails (stop()) if any REQUIRED row is FAIL.
#
# Value checks:
#   For each ome, reads qc_norm, feature_metadata, and sample_metadata from
#   staging_data_cache and runs compare_qc_norm() against the installed
#   MotrpacHumanPreSuspensionData package. Results print to console.
#
# Requires: config.R sourced; required_structure.R, expected_columns.R,
#           qc_norm_visualization_helpers.R sourced automatically below.

library(dplyr)
library(magrittr)
library(MotrpacHumanPreSuspensionData)

source(file.path(here::here(), "data-raw", "google_cloud_bucket_checks", "required_structure.R"))
source(file.path(here::here(), "data-raw", "google_cloud_bucket_checks", "expected_columns.R"))
source(file.path(here::here(), "data-raw", "google_cloud_bucket_checks", "qc_norm_visualization_helpers.R"))

# ---- List all files in staging -----------------------------------------------

message("Listing staging bucket: ", STAGING_BUCKET)
staging_files = system(paste(GSUTIL, "ls -R", STAGING_BUCKET), intern = TRUE)
staging_files = staging_files[!grepl(":$", staging_files)]
staging_files = staging_files[grepl("\\.txt$|\\.csv$|\\.txt\\.gz$|\\.html$", staging_files)]
message("Found ", length(staging_files), " files in staging.")

# ---- Shared cache for downloaded data frames ---------------------------------
# Keys: "<ome>__<tissue>__<type>" where type is qc_norm, feature_metadata, or sample_metadata.

staging_data_cache = new.env(parent = emptyenv())

# ---- Check QC report HTML presence per GCS subdir ---------------------------

.has_qc_report = function(subdir, staging_files) {
  pattern = paste0(STAGING_BUCKET, "/", subdir, "/.*\\.html$")
  any(grepl(pattern, staging_files))
}

# ---- Download a file to a temp location and return its path ------------------
# Caller owns cleanup (on.exit or explicit unlink).
# If you're doing a lot of testing, this can be substituted with a `dl_read_gcp`
# function to avoid re-downloading from the bucket each time.

.download_file = function(gcs_path) {
  tmp = tempfile(fileext = ".txt")
  system(paste(GSUTIL, "cp", shQuote(gcs_path), shQuote(tmp)), ignore.stdout = TRUE)
  if (!file.exists(tmp)) return(NULL)
  return(tmp)
}

# ---- Validate each manifest row ----------------------------------------------

validate_row = function(row) {
  ome           = row$ome
  tissue        = row$tissue
  tissue_code   = row$tissue_code
  if(ome == "clinical-chemistry") tissue_code = "t02-plasma"
  data_category = row$data_category
  data_details  = row$data_details
  gcs_subdir    = row$gcs_subdir

  # Find files that contain each required field (grep each independently, exact match).
  # Wrap data_details with underscores so "samples" does not match "removed-samples" paths.
  terms   = c(FILE_HEADER, tissue_code, ome, data_category, paste0("_", data_details, "_"))
  matched = staging_files
  for (term in terms) {
    matched = matched[grepl(term, matched, fixed = TRUE)]
  }

  if (length(matched) == 0) {
    if (data_details == "removed-samples") {
      return(data.frame(
        ome = ome, tissue = tissue, data_category = data_category,
        data_details = data_details, status = "PASS",
        reason = "No file — no removed samples",
        file_path = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    return(data.frame(
      ome = ome, tissue = tissue, data_category = data_category,
      data_details = data_details, status = "FAIL",
      reason = "No matching file found in staging",
      file_path = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  # Extract versions to select the best file
  ver_strings = regmatches(matched, regexpr("v(\\d+\\.\\d+)\\.txt", matched))
  versions    = as.numeric(sub("v", "", sub("\\.txt", "", ver_strings)))
  max_ver     = max(versions, na.rm = TRUE)
  best_file   = matched[which.max(versions)]

  # Determine whether this row's file should be cached for value comparison
  cache_type = if (data_category == "qc-norm") {
    "qc_norm"
  } else if (data_category == "metadata" && data_details == "features") {
    "feature_metadata"
  } else if (data_category == "metadata" && data_details == "samples") {
    "sample_metadata"
  } else {
    NULL
  }

  # Download the file once; on.exit ensures cleanup even on early return
  local_path = .download_file(best_file)
  if (is.null(local_path)) {
    return(data.frame(
      ome = ome, tissue = tissue, data_category = data_category,
      data_details = data_details, status = "FAIL",
      reason = "File found in staging but download failed",
      file_path = best_file,
      stringsAsFactors = FALSE
    ))
  }
  on.exit(unlink(local_path))

  # Column schema check
  expected = .get_expected_cols(data_category, ome, data_details)

  if (!is.null(expected) && length(expected$required) > 0) {
    actual_cols = tryCatch(
      colnames(read.csv(local_path, header = TRUE, sep = "\t", nrows = 1, check.names = FALSE)),
      error = function(e) character(0)
    )
    missing = setdiff(expected$required, actual_cols)
    extra   = intersect(expected$forbidden, actual_cols)

    if (length(missing) > 0 || length(extra) > 0) {
      reason_parts = character(0)
      if (length(missing) > 0) reason_parts = c(reason_parts, paste("missing cols:", paste(missing, collapse = ", ")))
      if (length(extra)   > 0) reason_parts = c(reason_parts, paste("forbidden cols present:", paste(extra, collapse = ", ")))
      return(data.frame(
        ome = ome, tissue = tissue, data_category = data_category,
        data_details = data_details, status = "FAIL",
        reason = paste(reason_parts, collapse = "; "),
        file_path = best_file,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Cache full data frame for value comparison (reuses the already-downloaded file)
  if (!is.null(cache_type)) {
    key = paste(ome, tissue, cache_type, sep = "__")
    staging_data_cache[[key]] = read_tsv(local_path)
  }

  return(data.frame(
    ome = ome, tissue = tissue, data_category = data_category,
    data_details = data_details, status = "PASS",
    reason = paste0("File found at v", max_ver),
    file_path = best_file,
    stringsAsFactors = FALSE
  ))
}

# ---- Run all rows ------------------------------------------------------------

message("Validating ", nrow(required_structure), " manifest entries...")
validation_results = do.call(rbind, lapply(seq_len(nrow(required_structure)), function(i) {
  validate_row(required_structure[i, ])
}))

# ---- Value checks: compare cached staging data vs installed package ----------

PROTEOMICS_OMES = c("prot-pr", "prot-ph", "prot-ol")

message("\n---- Value checks (staging vs installed package) ----")
existing_pkg = MotrpacHumanPreSuspensionData::load_qc(epigen = TRUE,
                                                      repo_local_dir = file.path(PRECOVID_REPO_PATH, "data", "tmp"),
                                                      remove_redundant_metab = FALSE,
                                                      load_acute_only = FALSE)

# Columns sourced from the pheno object are managed upstream; exclude them from
# the sample_metadata value comparison so they are never reported as differences.
pheno_cols = colnames(MotrpacHumanPreSuspensionData::load_pheno(load_acute_only = FALSE)[["pheno_data"]])

qc_ome_tissues = required_structure %>%
  dplyr::filter(data_category == "qc-norm") %>%
  dplyr::select(ome, tissue) %>%
  dplyr::distinct()

value_check_messages = character(0)
withCallingHandlers({
  for (ome_i in unique(qc_ome_tissues$ome)) {
    tissues_i = qc_ome_tissues$tissue[qc_ome_tissues$ome == ome_i]
    regen = list()
    for (tissue_i in tissues_i) {
      prefix = paste(ome_i, tissue_i, sep = "__")
      qc = staging_data_cache[[paste0(prefix, "__qc_norm")]]
      fm = staging_data_cache[[paste0(prefix, "__feature_metadata")]]
      sm = staging_data_cache[[paste0(prefix, "__sample_metadata")]]
      if (is.null(qc) || is.null(sm)) {
        message("  [", ome_i, "/", tissue_i, "] skipping value check: qc_norm or sample_metadata missing from cache")
        next
      }
      if (SKIP_TISSUE_NOT_IN_PKG && is.null(existing_pkg[[tissue_i]][[ome_i]])) {
        message("  [", ome_i, "/", tissue_i, "] skipping value check: loaded from cache but not present in installed package")
        next
      }
      if (is.null(fm)) {
        message("  [", ome_i, "/", tissue_i, "] note: no feature_metadata cached, feature_id check will be skipped")
      }
      regen[[tissue_i]] = list(qc_norm = qc, feature_metadata = fm, sample_metadata = sm)
    }
    if (length(regen) == 0) {
      message("[", ome_i, "] no tissues loaded from cache, skipping value check")
      next
    }
    fkey = if (ome_i %in% PROTEOMICS_OMES) "protein_id" else "feature_id"
    compare_qc_norm(ome_i, existing_pkg, regen, tissues = names(regen), feature_key = fkey,
                    pheno_cols = pheno_cols)
  }
}, message = function(m) {
  value_check_messages <<- c(value_check_messages, conditionMessage(m))
})

# ---- Check QC report presence per subdir ------------------------------------

missing_qc_subdirs = required_structure %>%
  dplyr::select(gcs_subdir) %>%
  dplyr::distinct() %>%
  dplyr::filter(!sapply(gcs_subdir, .has_qc_report, staging_files = staging_files)) %>%
  dplyr::pull(gcs_subdir)

qc_report_messages = if (length(missing_qc_subdirs) > 0) {
  paste0("WARN: No HTML QC report found in ", missing_qc_subdirs, "\n")
} else {
  character(0)
}

# ---- Write log ---------------------------------------------------------------

timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")
val_log   = file.path(LOG_DIR, paste0("structure_validation_", timestamp, ".tsv"))
write.table(validation_results, file = val_log, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n# QC reports\n", paste(qc_report_messages, collapse = ""),
    "\n# Value checks\n", paste(value_check_messages, collapse = ""),
    file = val_log, append = TRUE, sep = "")
message("Validation log written to: ", val_log)

# ---- Summary -----------------------------------------------------------------

status_counts = validation_results %>% dplyr::count(status)
message("Validation summary:")
print(status_counts)

required_fails = validation_results %>%
  dplyr::inner_join(
    required_structure %>% dplyr::filter(required) %>%
      dplyr::select(ome, tissue, data_category, data_details),
    by = c("ome", "tissue", "data_category", "data_details")
  ) %>%
  dplyr::filter(status == "FAIL")

if (nrow(required_fails) > 0) {
  message("\nFailed required entries:")
  print(required_fails)
  stop(nrow(required_fails), " required manifest entries FAILED. ",
       "Inspect validation log before proceeding: ", val_log)
}

message("Step 5 complete.")
