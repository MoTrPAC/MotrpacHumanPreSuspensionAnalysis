# Step 4: Apply the diff to the staging bucket.
# Uploads ADDED and MODIFIED files.
# For MODIFIED files, removes the old versioned file from staging after the new
# version is uploaded and its MD5 is verified.
# Verifies post-upload MD5 for every uploaded file before any deletion occurs.
# Writes logs/upload_log_{timestamp}.tsv.
#
# Requires: config.R sourced, Step 3 complete (diffs/latest.tsv)

library(dplyr)
library(magrittr)
library(tools)

latest_diff = file.path(DIFF_DIR, "latest.tsv")
if (!file.exists(latest_diff)) stop("diffs/latest.tsv not found. Run Step 3 first.")

diff_df = read.table(latest_diff, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

to_upload = diff_df %>% dplyr::filter(change_type %in% c("ADDED", "MODIFIED"))

message("Files to upload: ", nrow(to_upload))
message("  ADDED:    ", sum(to_upload$change_type == "ADDED"))
message("  MODIFIED: ", sum(to_upload$change_type == "MODIFIED"), " (old version will be removed after upload)")

upload_results = vector("list", nrow(to_upload))

# ---- Uploads -----------------------------------------------------------------

for (i in seq_len(nrow(to_upload))) {
  row = to_upload[i, ]
  exit = system(paste(GSUTIL, "cp", shQuote(row$local_path), shQuote(row$gcs_path)))

  # Post-upload MD5 verification
  stat_out  = system(paste(GSUTIL, "stat", shQuote(row$gcs_path)), intern = TRUE)
  md5_line  = stat_out[grepl("Hash \\(md5\\)", stat_out, ignore.case = TRUE)]
  remote_md5 = if (length(md5_line) > 0) trimws(sub(".*:\\s*", "", md5_line[1])) else NA_character_

  # gsutil returns base64 MD5; tools::md5sum returns hex — convert for comparison
  # gsutil stat md5 is base64-encoded; decode and convert to hex for comparison
  remote_md5_hex = tryCatch({
    raw_bytes = base64enc::base64decode(remote_md5)
    paste(sprintf("%02x", as.integer(raw_bytes)), collapse = "")
  }, error = function(e) NA_character_)

  verified = !is.na(remote_md5_hex) && (remote_md5_hex == row$new_md5)

  if (!verified) {
    message("MD5 MISMATCH: ", row$gcs_path)
    message("  local:  ", row$new_md5)
    message("  remote: ", remote_md5_hex)
  }

  # For MODIFIED files, remove the old versioned file from staging only — never
  # from the production snapshot. Derive the staging path by substituting the
  # production bucket prefix with STAGING_BUCKET. Only runs after verification.
  old_version_removed = NA
  if (row$change_type == "MODIFIED" && verified && !is.na(row$production_gcs_path)) {
    staging_old_path = sub(
      paste0("^", PRODUCTION_BUCKET, "/?"),
      paste0(STAGING_BUCKET, "/"),
      row$production_gcs_path
    )
    rm_exit = system(paste(GSUTIL, "rm", shQuote(staging_old_path)))
    old_version_removed = (rm_exit == 0)
    if (!old_version_removed) {
      message("WARNING: failed to remove old staging version: ", staging_old_path)
    }
  }

  upload_results[[i]] = data.frame(
    gcs_path            = row$gcs_path,
    production_gcs_path   = row$production_gcs_path,
    change_type         = row$change_type,
    old_version         = row$old_version,
    new_version         = row$new_version,
    exit_code           = exit,
    local_md5           = row$new_md5,
    remote_md5          = remote_md5_hex,
    verified            = verified,
    old_version_removed = old_version_removed,
    stringsAsFactors    = FALSE
  )
}

# ---- Write log ---------------------------------------------------------------

all_results = do.call(rbind, upload_results)
timestamp   = format(Sys.time(), "%Y%m%d_%H%M%S")
log_path    = file.path(LOG_DIR, paste0("upload_log_", timestamp, ".tsv"))
write.table(all_results, file = log_path, sep = "\t", row.names = FALSE, quote = FALSE)
message("Upload log written to: ", log_path)

upload_failures = all_results %>% dplyr::filter(!verified)
removal_failures = all_results %>%
  dplyr::filter(change_type == "MODIFIED", !is.na(old_version_removed), !old_version_removed)

if (nrow(upload_failures) > 0) {
  stop(nrow(upload_failures), " file(s) failed MD5 verification. ",
       "Inspect upload log before proceeding: ", log_path)
}

if (nrow(removal_failures) > 0) {
  stop(nrow(removal_failures), " old version(s) could not be removed from staging. ",
       "Inspect upload log before proceeding: ", log_path)
}

message("Step 4 complete: all uploads verified and old versions removed.")
