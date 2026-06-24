# Top-level driver for the GCS validation pipeline.
# Sources steps 01–08 in order. Run from the repo root or via here::here().
#
# Usage:
#   Rscript data-raw/google_cloud_bucket_checks/run_validation_pipeline.R
# or inside an R session:
#   source("data-raw/google_cloud_bucket_checks/run_validation_pipeline.R")

library(here)

pipeline_dir = file.path(here::here(), "data-raw", "google_cloud_bucket_checks")

step_files = c(
  "config.R",
  "01_snapshot_production.R",
  "02_copy_to_staging.R",
  "03_diff_local_vs_staging.R",
  "04_upload_to_staging.R",
  "05_validate_structure.R"
  #5 takes a long time
)

for (step in step_files) {
  step_path = file.path(pipeline_dir, step)
  message("\n", strrep("=", 60))
  message("Running: ", step)
  message(strrep("=", 60))
  source(step_path)
}

message("Pipeline complete. Review CHANGES.md and logs/ before merging branches.")
