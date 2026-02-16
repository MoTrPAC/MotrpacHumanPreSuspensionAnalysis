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
