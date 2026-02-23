test_that("run_cluster_cameraPR errors on invalid FCM input", {
  skip_if_not_installed("TMSig")
  expect_error(
    run_cluster_cameraPR(FCM = "not_a_list"),
    "output of run_cmeans"
  )
  expect_error(
    run_cluster_cameraPR(FCM = list(1, 2, 3)),
    "output of run_cmeans"
  )
})

test_that("run_cluster_ORA errors on invalid FCM input", {
  skip_if_not_installed("TMSig")
  expect_error(
    run_cluster_ORA(FCM = "not_a_list"),
    "output of run_cmeans"
  )
})
