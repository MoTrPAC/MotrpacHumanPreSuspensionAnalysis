# Pipeline configuration — edit these values before each release cycle.
# Sourced by every step script; do not run steps directly without sourcing this first.

# ---- Bucket paths ------------------------------------------------------------

PRODUCTION_BUCKET = "gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3"
STAGING_BUCKET    = "gs://pre-cawg/staging_20260511"

# ---- Version strings ---------------------------------------------------------

CURRENT_VERSION = "1.3"
NEW_VERSION     = "1.4"

# ---- Local paths -------------------------------------------------------------

 # Paths to downstream package repos (must match ~/config.json)
config_file = jsonlite::fromJSON("~/config.json")
PRECOVID_REPO_PATH     = config_file$precovid_repo_path
DATA_PKG_REPO_PATH     = config_file$data_package_path
ANALYSIS_PKG_REPO_PATH = here::here()

# Directory containing the incoming updated files to be diffed against staging.
# Subdirectory structure must mirror the GCS layout (epigenomics/, proteomics/, etc.)
LOCAL_UPDATED_DIR = file.path(PRECOVID_REPO_PATH,
                              "data", "tmp", "freeze")


# ---- Output paths (created automatically by each step) ----------------------

SNAPSHOT_DIR   = file.path(here::here(), "data-raw", "google_cloud_bucket_checks", "snapshots")
DIFF_DIR       = file.path(here::here(), "data-raw", "google_cloud_bucket_checks", "diffs")
LOG_DIR        = file.path(here::here(), "data-raw", "google_cloud_bucket_checks", "logs")
CHANGE_LOG_TSV = file.path(here::here(), "data-raw", "google_cloud_bucket_checks", "change_log.tsv")
CHANGES_MD     = file.path(here::here(), "data-raw", "google_cloud_bucket_checks", "CHANGES.md")

for (d in c(SNAPSHOT_DIR, DIFF_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- gsutil ------------------------------------------------------------------

GSUTIL = "gsutil"

# ---- BIC file name header ----------------------------------------------------

FILE_HEADER = "human-precovid-sed-adu"
