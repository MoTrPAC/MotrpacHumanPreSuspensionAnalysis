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

test_that("clinical chemistry is plotted only when a clinical ome is requested", {
  by_name <- plot_single_feature(
    feature = "Glucose",
    selected_omes = "metab-t-clinical",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_equal(unique(as.character(by_name$data$assay)), "metab-t-clinical")

  via_all <- plot_single_feature(
    feature = "Glucose",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_true("metab-t-clinical" %in% via_all$data$assay)

  # requesting another ome must not return clinical chemistry alongside it, and the
  # error has to name the ome that would have worked
  expect_error(plot_single_feature(feature = "Glucose",
                                   selected_omes = "transcript-rna-seq",
                                   selected_tissues = "blood",
                                   verbose = FALSE),
               "measured by clinical chemistry \\(metab-t-clinical\\)")

  # "metab" is the research platforms, and does not imply the clinical one. In blood the
  # analyte is also on the conventional panel, so the request is answerable without it —
  # what matters is that the clinical assay is not what comes back.
  research <- plot_single_feature(
    feature = "Glucose",
    selected_omes = "metab",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_false("metab-t-clinical" %in% research$data$assay)
  expect_true("metab-t-conv" %in% research$data$assay)

  # the clinical omes are separate from each other too
  expect_error(plot_single_feature(feature = "Glucose",
                                   selected_omes = "prot-clinical",
                                   selected_tissues = "blood",
                                   verbose = FALSE),
               "measured by clinical chemistry \\(metab-t-clinical\\)")

  # a feature that genuinely is not in the data keeps the original error
  expect_error(plot_single_feature(feature = "VEGFA",
                                   selected_omes = "metab-t-clinical",
                                   selected_tissues = "blood",
                                   verbose = FALSE),
               "No differential analysis corresponds")
})

test_that("the conventional metabolomics platform is plotted and labelled", {
  # metab-t-conv is no longer filtered out. It has no assay_codes row, so without a
  # fallback its facet strip reads NA.
  res <- plot_single_feature(
    feature = "Glucose",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_true("metab-t-conv" %in% res$data$assay)
  expect_false(any(is.na(res$data$tissue_assay)))
  expect_true("Blood Conv. Metab (log2)" %in% res$data$tissue_assay)
})

test_that("clinical prot analytes follow the same gate", {
  res <- plot_single_feature(
    feature = "Insulin",
    selected_omes = "prot-clinical",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_equal(unique(as.character(res$data$assay)), "prot-clinical")

  expect_error(plot_single_feature(feature = "Insulin",
                                   selected_omes = "prot-ol",
                                   selected_tissues = "blood",
                                   verbose = FALSE),
               "measured by clinical chemistry \\(prot-clinical\\)")
})

test_that("an analyte measured clinically and on a research platform splits by ome", {
  both <- plot_single_feature(
    feature = "Cortisol",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_true(all(c("metab-t-clinical", "metab-u-hilicpos") %in% both$data$assay))

  research_only <- plot_single_feature(
    feature = "Cortisol",
    selected_omes = "metab",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_false("metab-t-clinical" %in% research_only$data$assay)
  expect_true("metab-u-hilicpos" %in% research_only$data$assay)
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
