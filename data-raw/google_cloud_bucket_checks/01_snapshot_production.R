# Step 1: Snapshot the production bucket.
# Captures a timestamped MD5 inventory of every file in PRODUCTION_BUCKET.
# Writes snapshots/snapshot_{timestamp}.tsv and symlinks snapshots/latest.tsv.
#
# Requires: config.R sourced

library(dplyr)
library(magrittr)

timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")
snapshot_path = file.path(SNAPSHOT_DIR, paste0("snapshot_", timestamp, ".tsv"))
latest_path   = file.path(SNAPSHOT_DIR, "latest.tsv")

message("Listing all files in production bucket: ", PRODUCTION_BUCKET)
all_files = system(paste(GSUTIL, "ls -R", PRODUCTION_BUCKET), intern = TRUE)

# Keep only actual file paths (not directory listing headers ending in ":")
all_files = all_files[!grepl(":$", all_files)]
all_files = all_files[grepl("\\.txt$|\\.csv$|\\.txt\\.gz$|\\.rda$|\\.html$", all_files)]

if (length(all_files) == 0) stop("No files found in production bucket. Check PRODUCTION_BUCKET path and credentials.")

message("Found ", length(all_files), " files. Fetching MD5 and size via gsutil stat...")

snapshot_rows = lapply(all_files, function(gcs_path) {
  stat_out = system(paste(GSUTIL, "stat", shQuote(gcs_path)), intern = TRUE)

  md5_line  = stat_out[grepl("Hash \\(md5\\)", stat_out, ignore.case = TRUE)]
  size_line = stat_out[grepl("Content-Length", stat_out, ignore.case = TRUE)]

  md5  = if (length(md5_line)  > 0) trimws(sub(".*:\\s*", "", md5_line[1]))  else NA_character_
  size = if (length(size_line) > 0) trimws(sub(".*:\\s*", "", size_line[1])) else NA_character_

  data.frame(
    gcs_path           = gcs_path,
    size_bytes         = size,
    md5                = md5,
    snapshot_timestamp = timestamp,
    stringsAsFactors   = FALSE
  )
})

snapshot_df = do.call(rbind, snapshot_rows)

write.table(snapshot_df, file = snapshot_path, sep = "\t", row.names = FALSE, quote = FALSE)
message("Snapshot written to: ", snapshot_path)

# Overwrite the latest pointer
file.copy(snapshot_path, latest_path, overwrite = TRUE)
message("Latest snapshot pointer updated: ", latest_path)
message("Step 1 complete: ", nrow(snapshot_df), " files inventoried.")
