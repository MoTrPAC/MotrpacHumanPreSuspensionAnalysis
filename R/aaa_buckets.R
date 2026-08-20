# Bucket paths for the run-time epigenomics reads.
#
# precovid-repro is the source of truth: this mirrors STAGING_BUCKET in its
# config/pipeline.env. Bump both in the same release cycle, and keep it in step
# with the same constant in MotrpacHumanPreSuspensionData.
.STAGING_BUCKET <- "gs://pre-cawg/staging_20260806"
