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

test_that("plot_single_feature legends collect across tissues", {
  skip_if_not_installed("patchwork")

  tissues <- c("adipose", "blood", "muscle")
  plots <- lapply(tissues, function(tissue) {
    plot_single_feature(
      feature = "TAMALIN",
      selected_tissues = tissue,
      selected_omes = "transcript-rna-seq",
      verbose = FALSE
    )
  })

  # the reported failure needs a mix of tissues: at least one with no timepoint below
  # the p threshold, and at least one with some
  has_significant <- vapply(
    plots,
    function(p) any(p$data$below_p_cutoff == "Below p threshold"),
    logical(1)
  )
  skip_if_not(any(has_significant) && !all(has_significant),
              "TAMALIN no longer mixes significant and non-significant tissues")

  # every plot has to carry the same key set, or patchwork keeps one legend per variant
  fill_keys <- lapply(plots, ggplot2::get_guide_data, aesthetic = "fill")
  expect_equal(fill_keys[[2]], fill_keys[[1]])
  expect_equal(fill_keys[[3]], fill_keys[[1]])

  combined <- patchwork::wrap_plots(plots, ncol = 3, guides = "collect") &
    ggplot2::theme(legend.position = "right")
  assembled <- patchwork::patchworkGrob(combined)
  guide_box <- assembled$grobs[[which(assembled$layout$name == "guide-box")[1]]]

  # one legend for the exercise group colors, one for the p threshold fills
  expect_equal(sum(guide_box$layout$name == "guides"), 2)
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
