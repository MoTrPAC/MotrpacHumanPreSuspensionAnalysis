test_that("check_package_installation returns TRUE for installed package", {
  result <- MotrpacHumanPreSuspensionAnalysis:::check_package_installation(
    pkg = "stats", fun = "test_func"
  )
  expect_true(result)
})

test_that("check_package_installation errors for missing package", {
  expect_error(
    MotrpacHumanPreSuspensionAnalysis:::check_package_installation(
      pkg = "totally_fake_package_xyz123",
      fun = "some_function"
    ),
    "must be installed"
  )
})

test_that("check_package_installation includes task in error message", {
  expect_error(
    MotrpacHumanPreSuspensionAnalysis:::check_package_installation(
      pkg = "totally_fake_package_xyz123",
      fun = "some_function",
      task = "doing something"
    ),
    "doing something"
  )
})

test_that("check_package_installation gives special message for private data pkg", {
  skip_if(
    requireNamespace("MotrpacHumanPreSuspensionData", quietly = TRUE),
    "MotrpacHumanPreSuspensionData is installed; cannot test error path"
  )
  expect_error(
    MotrpacHumanPreSuspensionAnalysis:::check_package_installation(
      pkg = "MotrpacHumanPreSuspensionData",
      fun = "test"
    ),
    "private repository"
  )
})
