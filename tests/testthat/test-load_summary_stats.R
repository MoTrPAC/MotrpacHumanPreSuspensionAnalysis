test_that("load_summary_stats returns nested list by default", {
  res <- load_summary_stats(
    selected_omes = "prot-pr",
    selected_tissues = "muscle",
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_true("muscle" %in% names(res))
})

test_that("load_summary_stats single_matrix returns data.frame", {
  res <- load_summary_stats(
    selected_omes = "transcript-rna-seq",
    selected_tissues = "adipose",
    single_matrix = TRUE,
    verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
})

test_that("load_summary_stats handles 'all' tissues", {
  res <- load_summary_stats(
    selected_omes = "transcript-rna-seq",
    selected_tissues = "all",
    verbose = FALSE
  )
  tissues_returned <- names(res)
  expect_true(all(c("adipose", "blood", "muscle") %in% tissues_returned))
})

test_that("load_summary_stats rejects invalid tissue", {
  expect_error(
    load_summary_stats(selected_tissues = "liver", verbose = FALSE)
  )
})

test_that("load_summary_stats rejects invalid ome", {
  expect_error(
    load_summary_stats(selected_omes = "fake-ome", verbose = FALSE)
  )
})

test_that("metabolomics comes back as one stacked table per tissue", {
  res <- load_summary_stats(
    selected_omes = "metab",
    selected_tissues = "blood",
    verbose = FALSE
  )
  # One element named for the ome, not one per platform — the nesting
  # load_differential_analysis() returns.
  expect_named(res[["blood"]], "metab")

  stacked <- res[["blood"]][["metab"]]
  expect_true(all(stacked$assay == "metab"))
  expect_true("platform" %in% names(stacked))
  expect_gt(dplyr::n_distinct(stacked$platform), 1)
  # clinical chemistry is its own object and must not be in the stack
  expect_false("metab-t-clinical" %in% as.character(stacked$platform))
})

test_that("naming a single metabolomics platform returns the whole stack", {
  by_platform <- load_summary_stats(
    selected_omes = "metab-u-rppos",
    selected_tissues = "blood",
    verbose = FALSE
  )
  by_ome <- load_summary_stats(
    selected_omes = "metab",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_identical(by_platform, by_ome)
})

test_that("the stack is keyed the way its differential analysis is keyed", {
  sum_stats <- load_summary_stats(
    selected_omes = "metab", selected_tissues = "muscle", verbose = FALSE
  )[["muscle"]][["metab"]]
  da <- load_differential_analysis(
    selected_omes = "metab", selected_tissues = "muscle", verbose = FALSE
  )[["muscle"]][["metab"]]

  # The join that used to need a translation: both tiers now say assay = "metab" and name
  # the platform in `platform`, so (tissue, assay, platform, feature_id) means the same
  # thing on either side.
  expect_setequal(
    as.character(unique(sum_stats$platform)),
    as.character(unique(da$platform))
  )
  expect_true(all(unique(sum_stats$feature_id) %in% unique(as.character(da$feature_id))))
})

test_that("summary statistics and differential analysis nest the same way", {
  omes <- c("prot-pr", "metab")
  sum_stats <- load_summary_stats(selected_omes = omes, selected_tissues = "muscle",
                                  verbose = FALSE)
  da <- load_differential_analysis(selected_omes = omes, selected_tissues = "muscle",
                                   verbose = FALSE)
  expect_setequal(names(sum_stats[["muscle"]]), names(da[["muscle"]]))
})
