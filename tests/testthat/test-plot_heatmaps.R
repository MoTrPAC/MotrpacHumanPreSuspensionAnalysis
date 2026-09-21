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

# --- plot_feature_heatmap drawing options ---

heatmap_feature_ids <- c("ENSG00000109819.9", "ENSG00000112715.26",
                         "ENSG00000119508.18", "ENSG00000162772.17")

test_that("plot_feature_heatmap clusters rows across tissues", {
  skip_if_not_installed("TMSig")
  filename <- tempfile(fileext = ".pdf")
  expect_no_error(
    plot_feature_heatmap(
      feature_ids = heatmap_feature_ids,
      selected_tissue = c("muscle", "adipose"),
      selected_ome = "transcript-rna-seq",
      multi_tissue_clust_rows = TRUE,
      filename = filename
    )
  )
  expect_true(file.exists(filename))
})

test_that("plot_feature_heatmap return_drawing returns a drawing and page size", {
  skip_if_not_installed("TMSig")
  hm <- plot_feature_heatmap(
    feature_ids = heatmap_feature_ids,
    selected_tissue = "muscle",
    selected_ome = "transcript-rna-seq",
    return_drawing = TRUE
  )
  expect_named(hm, c("draw", "width", "height"))
  expect_type(hm$draw, "closure")
  expect_true(is.numeric(hm$width) && hm$width > 0)
  expect_true(is.numeric(hm$height) && hm$height > 0)

  filename <- tempfile(fileext = ".pdf")
  grDevices::pdf(filename, width = hm$width, height = hm$height)
  expect_no_error(hm$draw())
  grDevices::dev.off()
  expect_true(file.exists(filename))
})

test_that("plot_feature_heatmap passes row labels in matrix order to right_annotation", {
  skip_if_not_installed("TMSig")
  received <- NULL
  hm <- plot_feature_heatmap(
    feature_ids = heatmap_feature_ids,
    selected_tissue = c("muscle", "adipose"),
    selected_ome = "transcript-rna-seq",
    multi_tissue_clust_rows = TRUE,
    right_annotation = function(row_labels) {
      received <<- row_labels
      ComplexHeatmap::rowAnnotation(group = rep("A", length(row_labels)))
    },
    heatmap_args = list(row_names_side = "left"),
    return_drawing = TRUE
  )
  expect_setequal(received, c("PPARGC1A", "VEGFA", "NR4A3", "ATF3"))
  expect_identical(received, sort(received, method = "radix"))

  filename <- tempfile(fileext = ".pdf")
  grDevices::pdf(filename, width = hm$width, height = hm$height)
  expect_no_error(hm$draw())
  grDevices::dev.off()
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
