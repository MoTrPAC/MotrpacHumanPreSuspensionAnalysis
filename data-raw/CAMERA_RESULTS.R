CAMERA_RESULTS <- MotrpacHumanPreSuspensionAnalysis::run_cameraPR()

# Save
usethis::use_data(CAMERA_RESULTS,
                  internal = FALSE,
                  overwrite = TRUE,
                  compress = TRUE,
                  version = 3)
