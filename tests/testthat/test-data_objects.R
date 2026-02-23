# --- Color vectors ---

test_that("HUMAN_SEX_COLORS is a named character vector", {
  expect_type(HUMAN_SEX_COLORS, "character")
  expect_true(!is.null(names(HUMAN_SEX_COLORS)))
  expect_true(length(HUMAN_SEX_COLORS) >= 2L)
})

test_that("HUMAN_TISSUE_COLORS is a named character vector with expected names", {

  expect_type(HUMAN_TISSUE_COLORS, "character")
  expect_true(!is.null(names(HUMAN_TISSUE_COLORS)))
  expect_true(all(c("adipose", "blood", "muscle") %in% names(HUMAN_TISSUE_COLORS)))
})

test_that("HUMAN_TISSUE_ABBR is a named character vector", {
  expect_type(HUMAN_TISSUE_ABBR, "character")
  expect_true(!is.null(names(HUMAN_TISSUE_ABBR)))
  expect_true(length(HUMAN_TISSUE_ABBR) >= 3L)
})

test_that("HUMAN_EXERCISE_GROUP_COLORS is a named character vector", {
  expect_type(HUMAN_EXERCISE_GROUP_COLORS, "character")
  expect_true(!is.null(names(HUMAN_EXERCISE_GROUP_COLORS)))
  expect_true(length(HUMAN_EXERCISE_GROUP_COLORS) >= 2L)
})

test_that("HUMAN_ACUTE_TIMEPOINT_COLORS is a named character vector", {
  expect_type(HUMAN_ACUTE_TIMEPOINT_COLORS, "character")
  expect_true(!is.null(names(HUMAN_ACUTE_TIMEPOINT_COLORS)))
  expect_true(length(HUMAN_ACUTE_TIMEPOINT_COLORS) >= 3L)
})

test_that("HUMAN_OME_COLORS is a named character vector", {
  expect_type(HUMAN_OME_COLORS, "character")
  expect_true(!is.null(names(HUMAN_OME_COLORS)))
  expect_true(length(HUMAN_OME_COLORS) >= 4L)
})

# --- Reference tables ---

test_that("CONTRAST_CONVERTER has expected structure", {
  expect_s3_class(CONTRAST_CONVERTER, "data.frame")
  expect_equal(nrow(CONTRAST_CONVERTER), 33L)
  expected_cols <- c("contrast", "contrast_short", "contrast_type",
                     "contrast_category")
  expect_true(all(expected_cols %in% colnames(CONTRAST_CONVERTER)))
  expect_s3_class(CONTRAST_CONVERTER$contrast, "factor")
})

test_that("CONTRAST_CONVERTER contrast_type has exactly 5 levels", {
  types <- unique(as.character(CONTRAST_CONVERTER$contrast_type))
  expect_equal(
    sort(types),
    sort(c("exercise_with_controls", "exercise_no_controls",
           "Endur_vs_Resist", "baseline", "control_only"))
  )
})

test_that("HUMAN_FEATURE_TO_GENE has expected structure", {
  expect_s3_class(HUMAN_FEATURE_TO_GENE, "data.frame")
  expect_true(nrow(HUMAN_FEATURE_TO_GENE) > 1e6)
  expected_cols <- c("feature_id", "gene_symbol", "assay")
  expect_true(all(expected_cols %in% colnames(HUMAN_FEATURE_TO_GENE)))
})

test_that("SET_TO_ID has expected structure", {
  expect_s3_class(SET_TO_ID, "data.frame")
  expect_true(all(c("set_id", "set") %in% colnames(SET_TO_ID)))
  expect_true(nrow(SET_TO_ID) > 0)
})

test_that("OME_TISSUE_CODE exists and is a data.frame", {
  expect_s3_class(OME_TISSUE_CODE, "data.frame")
  expect_true(nrow(OME_TISSUE_CODE) > 0)
})

test_that("COVARIATES_FILE exists and is a data.frame", {
  expect_s3_class(COVARIATES_FILE, "data.frame")
  expect_true(nrow(COVARIATES_FILE) > 0)
})

test_that("OUTLIERS has expected structure", {
  expect_s3_class(OUTLIERS, "data.frame")
  expect_equal(nrow(OUTLIERS), 154L)
  expect_equal(ncol(OUTLIERS), 4L)
})

test_that("METABOLOMICS_CVS has expected structure", {
  expect_s3_class(METABOLOMICS_CVS, "data.frame")
  expect_true(nrow(METABOLOMICS_CVS) > 0)
  expected_cols <- c("feature_id", "tissue", "assay", "lowest_CV")
  expect_true(all(expected_cols %in% colnames(METABOLOMICS_CVS)))
})

test_that("UTORONTO_TFs exists and is a data.frame", {
  expect_s3_class(UTORONTO_TFs, "data.frame")
  expect_true(nrow(UTORONTO_TFs) > 0)
})

test_that("MOLECULAR_SIGNATURES is a named list of lists", {
  expect_type(MOLECULAR_SIGNATURES, "list")
  expect_true(!is.null(names(MOLECULAR_SIGNATURES)))
  expect_true(length(MOLECULAR_SIGNATURES) > 0)
  # Each top-level element should be a named list of character vectors
  expect_true(all(vapply(MOLECULAR_SIGNATURES, is.list, logical(1L))))
})

# --- Differential Analysis results ---

test_that("DA result objects exist and are data.frames with feature_id", {
  da_names <- c(
    "ADIPOSE_METAB_DA", "ADIPOSE_PROT_PH_DA", "ADIPOSE_PROT_PR_DA",
    "ADIPOSE_TRNSCRPT_DA", "BLOOD_METAB_DA", "BLOOD_PROT_OL_DA",
    "BLOOD_TRNSCRPT_DA", "MUSCLE_METAB_DA", "MUSCLE_PROT_PH_DA",
    "MUSCLE_PROT_PR_DA", "MUSCLE_TRNSCRPT_DA"
  )
  for (nm in da_names) {
    obj <- get(nm)
    expect_true(inherits(obj, "data.frame"), label = paste(nm, "is data.frame"))
    expect_true(nrow(obj) > 0, label = paste(nm, "has rows"))
    expect_true("feature_id" %in% colnames(obj),
                label = paste(nm, "has feature_id"))
  }
})

test_that("SPLICING_DA exists and has content", {
  expect_true(exists("SPLICING_DA"))
  expect_true(length(SPLICING_DA) > 0)
})

# --- Pre-computed enrichment and clustering results ---

test_that("CAMERA_RESULTS has expected structure", {
  expect_s3_class(CAMERA_RESULTS, "data.frame")
  expect_true(nrow(CAMERA_RESULTS) > 0)
  expected_cols <- c("tissue", "assay", "contrast", "set",
                     "p_value", "adj_p_value", "z.std")
  expect_true(all(expected_cols %in% colnames(CAMERA_RESULTS)))
})

test_that("FCM_CLUSTERS exists and has content", {
  expect_true(exists("FCM_CLUSTERS"))
  expect_true(length(FCM_CLUSTERS) > 0)
})

test_that("FCM_CAMERA has expected structure", {
  expect_s3_class(FCM_CAMERA, "data.frame")
  expect_true(nrow(FCM_CAMERA) > 0)
  expected_cols <- c("tissue", "assay", "cluster", "set", "p_value")
  expect_true(all(expected_cols %in% colnames(FCM_CAMERA)))
})

test_that("FCM_ORA has expected structure", {
  expect_s3_class(FCM_ORA, "data.frame")
  expect_true(nrow(FCM_ORA) > 0)
  expected_cols <- c("tissue", "assay", "cluster", "set", "p_value")
  expect_true(all(expected_cols %in% colnames(FCM_ORA)))
})
