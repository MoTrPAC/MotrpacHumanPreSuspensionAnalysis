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
  expect_length(res, 14L)
  expect_true(all(grepl("^metab-", res)))
  expect_false("transcript-rna-seq" %in% res)
  expect_false(any(duplicated(res)))
})

test_that("clinical_ome_list names the two v2.0 clinical omes", {
  res <- clinical_ome_list()
  expect_setequal(res, c("prot-clinical", "metab-t-clinical"))
  expect_true(all(res %in% ome_available_list()))
  # They are gated, not folded into the research platforms.
  expect_false("metab-t-clinical" %in% metab_only_list())
})

test_that("every loader takes load_clinical and defaults it to FALSE", {
  for (f in list(load_differential_analysis, load_summary_stats)) {
    expect_true("load_clinical" %in% names(formals(f)))
    expect_false(eval(formals(f)$load_clinical))
  }
})

test_that("load_clinical = FALSE keeps clinical chemistry out of 'all'", {
  # The regression this pins: adding the clinical omes to "all" silently changed
  # what every pre-v2.0 caller got. Four acute-repro panels failed on it — the
  # clinical DA rows overlap the combined metabolomics table on five analytes,
  # so a pivot keyed on (tissue, assay, feature_id) got duplicates.
  da <- suppressMessages(load_differential_analysis(
    selected_omes = "all", selected_tissues = "blood",
    single_matrix = TRUE, verbose = FALSE))
  expect_length(base::intersect(unique(da$assay), clinical_ome_list()), 0L)

  ss <- suppressMessages(load_summary_stats(
    selected_tissues = "blood", selected_omes = "all",
    single_matrix = TRUE, verbose = FALSE))
  expect_length(base::intersect(unique(ss$assay), clinical_ome_list()), 0L)
})

test_that("asking only for clinical with load_clinical = FALSE errors clearly", {
  # Without this the gate empties selected_omes and the failure surfaces much
  # later as a data.table error about a missing "contrast" column.
  expect_error(
    suppressMessages(load_differential_analysis(
      selected_omes = "prot-clinical", selected_tissues = "blood",
      single_matrix = TRUE, verbose = FALSE)),
    "load_clinical"
  )
  expect_error(
    suppressMessages(load_summary_stats(
      selected_tissues = "blood", selected_omes = "prot-clinical",
      single_matrix = TRUE, verbose = FALSE)),
    "load_clinical"
  )
})

test_that("load_clinical = FALSE drops clinical but keeps the rest", {
  da <- suppressMessages(load_differential_analysis(
    selected_omes = c("prot-ol", "prot-clinical"), selected_tissues = "blood",
    single_matrix = TRUE, verbose = FALSE))
  expect_equal(unique(da$assay), "prot-ol")
})

test_that("load_clinical = TRUE returns the clinical omes", {
  da <- suppressMessages(load_differential_analysis(
    selected_omes = "prot-clinical", selected_tissues = "blood",
    single_matrix = TRUE, verbose = FALSE, load_clinical = TRUE))
  expect_equal(unique(da$assay), "prot-clinical")

  ss <- suppressMessages(load_summary_stats(
    selected_tissues = "blood", selected_omes = "prot-clinical",
    single_matrix = TRUE, verbose = FALSE, load_clinical = TRUE))
  expect_equal(unique(ss$assay), "prot-clinical")
})

test_that("clinical metabolomics DA identifies itself rather than borrowing 'metab'", {
  # BLOOD_METAB_T_CLINICAL_DA stores assay = "metab", the same string the
  # combined table uses, which makes the two indistinguishable by
  # (tissue, assay, feature_id). The loader relabels on read.
  da <- suppressMessages(load_differential_analysis(
    selected_omes = "metab-t-clinical", selected_tissues = "blood",
    single_matrix = TRUE, verbose = FALSE, load_clinical = TRUE))
  expect_equal(unique(da$assay), "metab-t-clinical")
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
