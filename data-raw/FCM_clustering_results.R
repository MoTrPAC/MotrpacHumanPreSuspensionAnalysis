# Fuzzy c-means (FCM) clustering results
FCM_CLUSTERS <- MotrpacHumanPreSuspensionAnalysis::run_cmeans()

# CAMERA-PR applied to the cluster probabilities
FCM_CAMERA <- MotrpacHumanPreSuspensionAnalysis::run_cluster_cameraPR(
  FCM = FCM_CLUSTERS
)

# ORA applied to the disjoint clusters. These results are not used beyond being
# compared to the FCM_CAMERA results.
FCM_ORA <- MotrpacHumanPreSuspensionAnalysis::run_cluster_ORA(
  FCM = FCM_CLUSTERS
)

## Save ----
usethis::use_data(FCM_CLUSTERS,
                  internal = FALSE,
                  overwrite = TRUE,
                  compress = TRUE,
                  version = 3)

usethis::use_data(FCM_CAMERA,
                  internal = FALSE,
                  overwrite = TRUE,
                  compress = TRUE,
                  version = 3)

usethis::use_data(FCM_ORA,
                  internal = FALSE,
                  overwrite = TRUE,
                  compress = TRUE,
                  version = 3)
