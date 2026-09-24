test_that("run_cluster_cameraPR errors on invalid FCM input", {
  skip_if_not_installed("TMSig")
  expect_error(
    run_cluster_cameraPR(FCM = "not_a_list"),
    "named list of fclust objects"
  )
  expect_error(
    run_cluster_cameraPR(FCM = list(1, 2, 3)),
    "named list of fclust objects"
  )
})

test_that("run_cluster_ORA errors on invalid FCM input", {
  skip_if_not_installed("TMSig")
  expect_error(
    run_cluster_ORA(FCM = "not_a_list"),
    "named list of fclust objects"
  )
})
