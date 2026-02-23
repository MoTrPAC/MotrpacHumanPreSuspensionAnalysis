test_that("ome_available_list returns correct character vector", {
  res <- ome_available_list()
  expect_type(res, "character")
  expect_length(res, 21L)
  expect_true("transcript-rna-seq" %in% res)
  expect_true("prot-pr" %in% res)
  expect_true("prot-ph" %in% res)
  expect_true("prot-ol" %in% res)
  expect_true("metab-u-hilicpos" %in% res)
  expect_true("epigen-atac-seq" %in% res)
  expect_true("epigen-methylcap-seq" %in% res)
  expect_false(any(duplicated(res)))
})

test_that("metab_only_list returns correct character vector", {
  res <- metab_only_list()
  expect_type(res, "character")
  expect_length(res, 14L)
  expect_true(all(grepl("^metab-", res)))
  expect_false("transcript-rna-seq" %in% res)
  expect_false(any(duplicated(res)))
})

test_that("metab_only_list is a subset of ome_available_list", {
  expect_true(all(metab_only_list() %in% ome_available_list()))
})

test_that("tissue_available_list returns correct character vector", {
  res <- suppressMessages(tissue_available_list())
  expect_type(res, "character")
  expect_length(res, 3L)
  expect_equal(sort(res), c("adipose", "blood", "muscle"))
})

test_that("tissue_available_list emits message when verbose = TRUE", {
  expect_message(tissue_available_list(verbose = TRUE))
})

test_that("tissue_available_list suppresses message when verbose = FALSE", {
  expect_no_message(tissue_available_list(verbose = FALSE))
})

# --- Internal helpers ---

test_that(".find_ome finds ome in file path", {
  fn <- MotrpacHumanPreSuspensionAnalysis:::.find_ome
  expect_equal(
    fn("gs://bucket/t02-transcript-rna-seq/file.txt"),
    "transcript-rna-seq"
  )
  expect_equal(fn("gs://bucket/t06-prot-pr/file.txt"), "prot-pr")
  expect_equal(fn("something/metab-u-hilicpos/data.csv"), "metab-u-hilicpos")
})

test_that(".find_ome returns NULL when no ome matches", {
  fn <- MotrpacHumanPreSuspensionAnalysis:::.find_ome
  expect_null(fn("gs://bucket/no_ome_here/file.txt"))
})

test_that(".find_tissue finds tissue from code", {
  fn <- MotrpacHumanPreSuspensionAnalysis:::.find_tissue
  expect_equal(fn("gs://bucket/t02-something"), "blood")
  expect_equal(fn("gs://bucket/t06-something"), "muscle")
  expect_equal(fn("gs://bucket/t07-something"), "adipose")
  expect_equal(fn("gs://bucket/t10-something"), "muscle")
  expect_equal(fn("gs://bucket/t11-something"), "adipose")
})

test_that(".find_tissue returns NULL for no match", {
  fn <- MotrpacHumanPreSuspensionAnalysis:::.find_tissue
  expect_null(fn("gs://bucket/t99-something"))
})
