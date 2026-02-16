test_that("load_differential_analysis returns nested list by default", {
  res <- load_differential_analysis(
    selected_omes = "prot-pr",
    selected_tissues = "muscle",
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_true("muscle" %in% names(res))
  expect_true("prot-pr" %in% names(res[["muscle"]]))
  expect_s3_class(res[["muscle"]][["prot-pr"]], "data.table")
})

test_that("load_differential_analysis single_matrix returns data.table", {
  res <- load_differential_analysis(
    selected_omes = "prot-pr",
    selected_tissues = "muscle",
    single_matrix = TRUE,
    verbose = FALSE
  )
  expect_s3_class(res, "data.table")
  expect_true(nrow(res) > 0)
  expected_cols <- c("tissue", "assay", "feature_id", "logFC",
                     "contrast", "p_value", "adj_p_value")
  expect_true(all(expected_cols %in% colnames(res)))
})

test_that("load_differential_analysis handles 'all' tissues", {
  res <- load_differential_analysis(
    selected_omes = "transcript-rna-seq",
    selected_tissues = "all",
    verbose = FALSE
  )
  expect_true(all(c("adipose", "blood", "muscle") %in% names(res)))
})

test_that("load_differential_analysis handles multiple omes", {
  res <- load_differential_analysis(
    selected_omes = c("prot-pr", "transcript-rna-seq"),
    selected_tissues = "muscle",
    verbose = FALSE
  )
  expect_true("muscle" %in% names(res))
  expect_true(all(c("prot-pr", "transcript-rna-seq") %in%
                     names(res[["muscle"]])))
})

test_that("load_differential_analysis rejects invalid tissue", {
  expect_error(
    load_differential_analysis(selected_tissues = "liver", verbose = FALSE)
  )
})

test_that("load_differential_analysis rejects invalid ome", {
  expect_error(
    load_differential_analysis(selected_omes = "fake-ome", verbose = FALSE)
  )
})

test_that("load_differential_analysis metab consolidates to single key", {
  res <- load_differential_analysis(
    selected_omes = "metab-u-hilicpos",
    selected_tissues = "muscle",
    verbose = FALSE
  )
  expect_true("muscle" %in% names(res))
  expect_true("metab" %in% names(res[["muscle"]]))
})

test_that("load_differential_analysis combine_with_featgene adds gene columns", {
  res <- load_differential_analysis(
    selected_omes = "transcript-rna-seq",
    selected_tissues = "blood",
    combine_with_featgene = TRUE,
    verbose = FALSE
  )
  dt <- res[["blood"]][["transcript-rna-seq"]]
  expect_true("gene_symbol" %in% colnames(dt))
})

test_that("load_differential_analysis DA tables have expected column types", {
  res <- load_differential_analysis(
    selected_omes = "prot-pr",
    selected_tissues = "muscle",
    single_matrix = TRUE,
    verbose = FALSE
  )
  expect_type(res$logFC, "double")
  expect_type(res$p_value, "double")
  expect_type(res$adj_p_value, "double")
  expect_s3_class(res$contrast, "factor")
})
