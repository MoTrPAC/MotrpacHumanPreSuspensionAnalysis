test_that("run_cmeans returns named list with fclust objects", {
  skip_if_not_installed("Mfuzz")
  skip_if_not_installed("Biobase")

  res <- run_cmeans(
    selected_tissues = "blood",
    selected_omes = "transcript-rna-seq",
    num_clusters_blood = 3L
  )
  expect_type(res, "list")
  expect_true("blood" %in% names(res))

  fclust_obj <- res[["blood"]]
  expect_s3_class(fclust_obj, "fclust")
  expect_true("membership" %in% names(fclust_obj))
  expect_true("centers" %in% names(fclust_obj))
  expect_true("cluster" %in% names(fclust_obj))
  expect_true("input" %in% names(fclust_obj))
  expect_equal(nrow(fclust_obj[["centers"]]), 3L)
})
