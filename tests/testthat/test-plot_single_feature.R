test_that("plot_single_feature returns ggplot for a known gene", {
  res <- plot_single_feature(
    feature = "VEGFA",
    selected_tissues = "muscle",
    selected_omes = "transcript-rna-seq",
    verbose = FALSE
  )
  expect_s3_class(res, "gg")
  expect_s3_class(res, "ggplot")
})

test_that("plot_single_feature returns ggplot for a known metabolite", {
  res <- plot_single_feature(
    feature = "CAR 10:0",
    selected_omes = "metab",
    verbose = FALSE
  )
  expect_s3_class(res, "gg")
  expect_s3_class(res, "ggplot")
})

test_that("plot_single_feature returns ggplot for a clinical analyte", {
  res <- plot_single_feature(
    feature = "Glucose",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_s3_class(res, "gg")
  expect_s3_class(res, "ggplot")
})

test_that("plot_single_feature rejects invalid tissue", {
  expect_error(
    plot_single_feature(
      feature = "VEGFA",
      selected_tissues = "liver"
    )
  )
})

test_that("plot_single_feature rejects invalid ome", {
  expect_error(
    plot_single_feature(
      feature = "VEGFA",
      selected_omes = "fake-ome"
    )
  )
})

test_that("plot_single_feature saves to file when output_file provided", {
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)

  res <- plot_single_feature(
    feature = "VEGFA",
    selected_tissues = "muscle",
    selected_omes = "transcript-rna-seq",
    output_file = tmp,
    verbose = FALSE
  )
  expect_true(file.exists(tmp))
  expect_s3_class(res, "ggplot")
})
