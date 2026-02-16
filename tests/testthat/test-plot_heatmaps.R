# --- plot_feature_heatmap input validation ---

test_that("plot_feature_heatmap validates post_min correctly", {
  expect_error(
    plot_feature_heatmap(
      feature_ids = "TP53",
      selected_tissue = "muscle",
      selected_ome = "prot-pr",
      post_min = 99,
      filename = tempfile(fileext = ".pdf")
    ),
    "15, 30, or 45"
  )
})

test_that("plot_feature_heatmap validates post_min is numeric", {
  expect_error(
    plot_feature_heatmap(
      feature_ids = "TP53",
      selected_tissue = "muscle",
      selected_ome = "prot-pr",
      post_min = "abc",
      filename = tempfile(fileext = ".pdf")
    ),
    "numeric"
  )
})

test_that("plot_feature_heatmap validates post_hr correctly", {
  expect_error(
    plot_feature_heatmap(
      feature_ids = "TP53",
      selected_tissue = "muscle",
      selected_ome = "prot-pr",
      post_hr = 5,
      filename = tempfile(fileext = ".pdf")
    ),
    "3.5 or 4"
  )
})

test_that("plot_feature_heatmap rejects multiple omes", {
  expect_error(
    plot_feature_heatmap(
      feature_ids = "TP53",
      selected_tissue = "muscle",
      selected_ome = c("prot-pr", "prot-ph"),
      filename = tempfile(fileext = ".pdf")
    ),
    "only 1 ome"
  )
})

test_that("plot_feature_heatmap validates full_modality_names is logical", {
  expect_error(
    plot_feature_heatmap(
      feature_ids = "TP53",
      selected_tissue = "muscle",
      selected_ome = "prot-pr",
      full_modality_names = "yes",
      filename = tempfile(fileext = ".pdf")
    ),
    "logical"
  )
})

# --- plot_enrich_heatmap input validation ---

test_that("plot_enrich_heatmap errors on missing required columns", {
  skip_if_not_installed("TMSig")
  expect_error(
    plot_enrich_heatmap(
      x = data.frame(a = 1),
      selected_ome = "prot-pr",
      filename = tempfile(fileext = ".pdf")
    ),
    "missing the following required column"
  )
})

test_that("plot_enrich_heatmap errors when x lacks z.std or NES", {
  skip_if_not_installed("TMSig")
  fake_df <- data.frame(
    tissue = "muscle",
    assay = "prot-pr",
    contrast = "c1",
    contrast_type = "exercise_with_controls",
    set_id = "00001",
    set = "s1",
    set_short = "s1",
    p_value = 0.01,
    adj_p_value = 0.05
  )
  expect_error(
    plot_enrich_heatmap(
      x = fake_df,
      selected_ome = "prot-pr",
      filename = tempfile(fileext = ".pdf")
    ),
    "z.std.*NES"
  )
})

# --- plot_cluster_enrichment input validation ---

test_that("plot_cluster_enrichment requires PDF filename", {
  skip_if_not_installed("TMSig")
  expect_error(
    plot_cluster_enrichment(
      x = FCM_CAMERA,
      filename = "not_a_pdf.txt"
    ),
    "PDF"
  )
})

test_that("plot_cluster_enrichment errors on missing columns", {
  skip_if_not_installed("TMSig")
  expect_error(
    plot_cluster_enrichment(
      x = data.frame(a = 1),
      filename = tempfile(fileext = ".pdf")
    ),
    "missing the following required column"
  )
})
