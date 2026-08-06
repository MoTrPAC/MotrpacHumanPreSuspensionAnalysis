test_that("ome_available_list returns correct character vector", {
  res <- ome_available_list()
  expect_type(res, "character")
  expect_length(res, 22L)
  expect_true("transcript-rna-seq" %in% res)
  expect_true("prot-pr" %in% res)
  expect_true("prot-ph" %in% res)
  expect_true("prot-ol" %in% res)
  expect_true("metab-u-hilicpos" %in% res)
  expect_true("epigen-atac-seq" %in% res)
  expect_true("epigen-methylcap-seq" %in% res)
  expect_false(any(duplicated(res)))
})

test_that("ome_available_list carries the v2.0 clinical split", {
  res <- ome_available_list()
  expect_true("metab-t-clinical" %in% res)
  expect_true("prot-clinical" %in% res)
  # v1.3 named clinical chemistry as one assay; v2.0 splits it in two.
  expect_false("clinical-chemistry" %in% res)
})

test_that("ome_available_list offers no ome the pipeline stopped building", {
  # metab-meta-reg was dropped from the pipeline entirely, so it produced no
  # object while every accessor still offered it as a choice that matched
  # nothing.
  expect_false("metab-meta-reg" %in% ome_available_list())
})

test_that("every ome this package ships is reachable through ome_available_list", {
  # The failure this pins is an object that ships but no accessor can return,
  # because its ome is absent from the vocabulary the loaders gate on. Omes are
  # derived here exactly the way load_summary_stats() and
  # load_differential_analysis() derive them, so the test moves if they do.
  items <- data(package = "MotrpacHumanPreSuspensionAnalysis")[["results"]][, "Item"]

  ome_of <- function(x, suffix) {
    inner <- sub(paste0("^[^_]+_(.*)_", suffix, "$"), "\\1", x)
    out <- gsub("_", "-", tolower(inner))
    out[out == "trnscrpt"] <- "transcript-rna-seq"
    out
  }

  shipped <- c(ome_of(grep("_SUM_STATS$", items, value = TRUE), "SUM_STATS"),
               ome_of(grep("_DA$", items, value = TRUE), "DA"))
  # Two names are not omes and are expected here:
  #   "metab"       the combined per-tissue DA table, not a platform
  #   "splicing-da" SPLICING_DA carries no tissue prefix, so the TISSUE_OME_DA
  #                 pattern leaves it unchanged. The loader derives its tissue as
  #                 "splicing", which matches no selectable tissue, so the object
  #                 is reached by name rather than through load_differential_analysis().
  shipped <- setdiff(unique(shipped), c("metab", "splicing-da"))

  expect_equal(sort(setdiff(shipped, ome_available_list())), character(0))
})

test_that("metab_only_list returns correct character vector", {
  res <- metab_only_list()
  expect_type(res, "character")
  expect_length(res, 15L)
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
