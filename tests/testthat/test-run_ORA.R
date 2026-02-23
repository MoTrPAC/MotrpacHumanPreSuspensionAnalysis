test_that("run_ORA returns data.frame with expected columns", {
  skip_if_not_installed("TMSig")

  bg <- MOLECULAR_SIGNATURES
  bg[c("PSP", "REFMET")] <- NULL
  bg <- unique(unlist(bg))

  input <- MOLECULAR_SIGNATURES$MITOCARTA[["MITOCARTA_OXPHOS"]]

  res <- run_ORA(
    input = input,
    background = bg,
    database = "REACTOME",
    min_size = 5L
  )
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
  expected_cols <- c("set", "set_size", "p_value", "adj_p_value",
                     "set_size_in_input", "input_size", "background_size")
  expect_true(all(expected_cols %in% colnames(res)))
})

test_that("run_ORA p-values are in valid range", {
  skip_if_not_installed("TMSig")

  bg <- MOLECULAR_SIGNATURES
  bg[c("PSP", "REFMET")] <- NULL
  bg <- unique(unlist(bg))

  input <- MOLECULAR_SIGNATURES$MITOCARTA[["MITOCARTA_OXPHOS"]]

  res <- run_ORA(input = input, background = bg, database = "REACTOME")
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
  expect_true(all(res$adj_p_value >= 0 & res$adj_p_value <= 1))
})

test_that("run_ORA errors on non-character input", {
  skip_if_not_installed("TMSig")
  expect_error(
    run_ORA(input = 1:5, background = letters),
    "character vector"
  )
})

test_that("run_ORA errors on empty input", {
  skip_if_not_installed("TMSig")
  expect_error(
    run_ORA(input = character(0), background = letters),
    "character vector"
  )
})

test_that("run_ORA errors on non-character background", {
  skip_if_not_installed("TMSig")
  expect_error(
    run_ORA(input = "A", background = 1:5),
    "character vector"
  )
})

test_that("run_ORA errors when input is not subset of background", {
  skip_if_not_installed("TMSig")
  expect_error(
    run_ORA(
      input = c("NOT_A_GENE_XYZ"),
      background = c("BRCA1", "TP53")
    ),
    "subset of"
  )
})
