# Step 2: Copy production bucket to staging folder.
# Uses gsutil rsync -r (single-threaded) to mirror PRODUCTION_BUCKET into STAGING_BUCKET,
# preserving the full subdirectory structure (epigenomics/, proteomics/, etc.).
# rsync is preferred over cp -r here because it preserves the folder hierarchy
# correctly when copying between GCS paths — cp with ** glob does not.
# Verifies that the file count in staging matches the snapshot from Step 1.
#
# Requires: config.R sourced, Step 1 complete (snapshots/latest.tsv)

library(dplyr)
library(magrittr)

latest_path = file.path(SNAPSHOT_DIR, "latest.tsv")
if (!file.exists(latest_path)) stop("snapshots/latest.tsv not found. Run Step 1 first.")

snapshot_df = read.table(latest_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
expected_count = nrow(snapshot_df)

message("Copying ", PRODUCTION_BUCKET, " -> ", STAGING_BUCKET)
message("This may take several minutes for large buckets.")

# NOTE: Do not add -m (parallel/multi-threading) flag here. It causes intermittent
# copy failures when rsync-ing between GCS buckets, likely due to rate limiting or
# connection contention on large transfers. Run single-threaded despite being slower.
exit_code = system(paste(
  GSUTIL, "rsync -r",
  shQuote(paste0(PRODUCTION_BUCKET, "/")),
  shQuote(paste0(STAGING_BUCKET, "/"))
))

if (exit_code != 0) stop("gsutil rsync failed with exit code: ", exit_code)

# Verify count in staging
message("Verifying file count in staging...")
staging_files = system(paste(GSUTIL, "ls -R", STAGING_BUCKET), intern = TRUE)
staging_files = staging_files[!grepl(":$", staging_files)]
staging_files = staging_files[grepl("\\.txt$|\\.csv$|\\.txt\\.gz$|\\.rda$|\\.html$", staging_files)]

actual_count = length(staging_files)
verified = (actual_count >= expected_count)

status_line = paste0("staging_copy_verified: ", verified,
                     " | expected: ", expected_count,
                     " | found: ", actual_count)
message(status_line)

log_path = file.path(LOG_DIR, paste0("copy_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
writeLines(c(status_line, paste("timestamp:", Sys.time())), log_path)

if (!verified) {
  stop("Staging file count (", actual_count, ") is less than expected (",
       expected_count, "). Inspect staging bucket before proceeding.")
}

message("Step 2 complete: staging bucket verified with ", actual_count, " files.")
