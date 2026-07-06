# Step 3: Diff local updated files against the staging snapshot.
# Compares LOCAL_UPDATED_DIR against snapshots/latest.tsv (production baseline).
#
# Two-phase exclusion logic:
#   Phase 1 — Files whose exact basename (including version) already exists in
#             staging are skipped entirely before any processing. These are
#             definitively unchanged and never enter the diff.
#   Phase 2 — Remaining files are joined by a version-stripped basename (join_key),
#             so a file bumped from v1.3 to v1.4 is matched to its staging
#             counterpart rather than appearing as ADDED + REMOVED.
#
# Change categories in output:
#   ADDED    — local file has no matching join_key in staging
#   MODIFIED — join_key matches; versions or content differ
#   (REMOVED entries — staging files with no local match — are excluded from output)
#
# Output columns in diff_df:
#   local_path, gcs_path, production_gcs_path, join_key,
#   ome, tissue_code, data_category, data_details,
#   old_version, new_version, old_md5, new_md5, change_type
#
# Writes diffs/diff_{timestamp}.tsv and diffs/latest.tsv.
# Requires: config.R sourced, Step 1 complete

library(dplyr)
library(magrittr)
library(tools)

latest_path = file.path(SNAPSHOT_DIR, "latest.tsv")
if (!file.exists(latest_path)) stop("snapshots/latest.tsv not found. Run Step 1 first.")

snapshot_df = read.table(latest_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (!dir.exists(LOCAL_UPDATED_DIR)) {
  stop("LOCAL_UPDATED_DIR does not exist: ", LOCAL_UPDATED_DIR)
}

# ---- List all local files ----------------------------------------------------

local_files = list.files(LOCAL_UPDATED_DIR, recursive = TRUE, full.names = TRUE)
local_files = local_files[grepl("\\.txt$|\\.csv$|\\.txt\\.gz$|\\.html$", local_files)]
message("Found ", length(local_files), " local files in ", LOCAL_UPDATED_DIR)

# ---- Exclude files with the same basename (same version) as staging ----------
# Basename encodes the version string, so an identical basename means no change.
# These are excluded from both sides of the join to avoid false REMOVED flags.

staging_basenames  = basename(snapshot_df$gcs_path)
same_basename      = basename(local_files) %in% staging_basenames
skipped_basenames  = basename(local_files)[same_basename]

if (any(same_basename)) {
  message("Skipping ", sum(same_basename), " file(s) whose basename already exists in staging (same version — no diff needed).")
  local_files = local_files[!same_basename]
}

# ---- Helpers -----------------------------------------------------------------

.local_to_gcs = function(local_path, local_base, staging_bucket) {
  rel = sub(paste0("^", normalizePath(local_base), "/?"), "", normalizePath(local_path))
  paste0(staging_bucket, "/", rel)
}

# File name pattern: human-precovid-sed-adu_{tissue_code}_{ome}_{category}_{details}_v{ver}.txt
.parse_file_name = function(fname) {
  base  = basename(fname)
  parts = strsplit(base, "_")[[1]]
  if (length(parts) < 6) return(list(ome = NA, tissue_code = NA, data_category = NA, data_details = NA, version = NA))

  tissue_code   = parts[2]
  ome           = parts[3]
  data_category = parts[4]
  data_details  = paste(parts[5:(length(parts) - 1)], collapse = "_")
  version_raw   = parts[length(parts)]
  version       = sub("^v", "", sub("\\.txt(\\.gz)?$", "", version_raw))

  return(list(
    ome           = ome,
    tissue_code   = tissue_code,
    data_category = data_category,
    data_details  = data_details,
    version       = version
  ))
}

# Strip the version suffix to produce a join key that matches across versions.
# e.g. human-precovid-sed-adu_t06-muscle_transcript-rna-seq_qc-norm_log-cpm_v1.3.txt
#   -> human-precovid-sed-adu_t06-muscle_transcript-rna-seq_qc-norm_log-cpm
.strip_version = function(fname) {
  sub("_v[0-9]+\\.[0-9]+\\.txt(\\.gz)?$", "", basename(fname))
}

# ---- Build local file table --------------------------------------------------

local_rows = lapply(local_files, function(lp) {
  meta = .parse_file_name(lp)
  data.frame(
    local_path    = lp,
    gcs_path      = .local_to_gcs(lp, LOCAL_UPDATED_DIR, STAGING_BUCKET),
    join_key      = .strip_version(lp),
    ome           = meta$ome,
    tissue_code   = meta$tissue_code,
    data_category = meta$data_category,
    data_details  = meta$data_details,
    new_version   = meta$version,
    new_md5       = as.character(md5sum(lp)),
    stringsAsFactors = FALSE
  )
})

local_df = do.call(rbind, local_rows)

# ---- Build snapshot table ----------------------------------------------------

snapshot_small = snapshot_df %>%
  dplyr::filter(!basename(gcs_path) %in% skipped_basenames) %>%
  dplyr::mutate(
    join_key    = .strip_version(gcs_path),
    old_version = sub(".*_v([0-9]+\\.[0-9]+)\\.txt(\\.gz)?$", "\\1", basename(gcs_path))
  ) %>%
  dplyr::select(join_key, production_gcs_path = gcs_path, old_version, old_md5 = md5)

# ---- Join on version-stripped basename and classify changes ------------------

diff_df = dplyr::full_join(local_df, snapshot_small, by = "join_key") %>%
  dplyr::mutate(
    change_type = dplyr::case_when(
      !is.na(new_md5) & is.na(old_md5)  ~ "ADDED",
      !is.na(new_md5) & !is.na(old_md5) ~ "MODIFIED",
      is.na(new_md5)  & !is.na(old_md5) ~ "REMOVED",
      TRUE                               ~ "UNCHANGED"
    )
  ) %>%
  dplyr::select(
    local_path, gcs_path, production_gcs_path, join_key,
    ome, tissue_code, data_category, data_details,
    old_version, new_version, old_md5, new_md5, change_type
  ) %>%
  dplyr::filter(!is.na(local_path))

# ---- Summarize and write -----------------------------------------------------

summary_tbl = diff_df %>%
  dplyr::count(change_type) %>%
  dplyr::arrange(change_type)

message("Diff summary:")
print(summary_tbl)

timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")
diff_path  = file.path(DIFF_DIR, paste0("diff_", timestamp, ".tsv"))
latest_diff = file.path(DIFF_DIR, "latest.tsv")

write.table(diff_df, file = diff_path, sep = "\t", row.names = FALSE, quote = FALSE)
file.copy(diff_path, latest_diff, overwrite = TRUE)

message("Diff written to: ", diff_path)
message("Step 3 complete.")
