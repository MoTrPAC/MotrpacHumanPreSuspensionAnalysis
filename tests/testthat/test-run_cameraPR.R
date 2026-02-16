test_that("run_cameraPR returns data.frame with expected columns", {
  skip_if_not_installed("TMSig")

  res <- run_cameraPR(
    selected_omes = "prot-pr",
    selected_tissues = "muscle",
    database = "REACTOME"
  )
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
  expected_cols <- c("tissue", "assay", "contrast_type", "contrast",
                     "set", "set_size", "direction", "z.std",
                     "p_value", "adj_p_value")
  expect_true(all(expected_cols %in% colnames(res)))
})

test_that("run_cameraPR direction is Up or Down", {
  skip_if_not_installed("TMSig")

  res <- run_cameraPR(
    selected_omes = "prot-pr",
    selected_tissues = "muscle",
    database = "REACTOME"
  )
  expect_true(all(as.character(res$direction) %in% c("Up", "Down")))
})

test_that("run_cameraPR p-values are in valid range", {
  skip_if_not_installed("TMSig")

  res <- run_cameraPR(
    selected_omes = "prot-pr",
    selected_tissues = "muscle",
    database = "REACTOME"
  )
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
  expect_true(all(res$adj_p_value >= 0 & res$adj_p_value <= 1))
})
