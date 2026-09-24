test_that("plot_single_feature returns ggplot for a known gene", {
  res <- plot_single_feature(
    feature = "VEGFA",
    selected_tissues = "muscle",
    selected_omes = "transcript-rna-seq",
    verbose = FALSE
  )
  expect_s3_class(res, "gg")
  expect_s3_class(res, "ggplot")
})

test_that("plot_single_feature returns ggplot for a known metabolite", {
  res <- plot_single_feature(
    feature = "CAR 10:0",
    selected_omes = "metab",
    verbose = FALSE
  )
  expect_s3_class(res, "gg")
  expect_s3_class(res, "ggplot")
})

test_that("plot_single_feature returns ggplot for a clinical analyte", {
  res <- plot_single_feature(
    feature = "Glucose",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_s3_class(res, "gg")
  expect_s3_class(res, "ggplot")
})

test_that("clinical chemistry is plotted only when a clinical ome is requested", {
  by_name <- plot_single_feature(
    feature = "Glucose",
    selected_omes = "metab-t-clinical",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_equal(unique(as.character(by_name$data$assay)), "metab-t-clinical")

  via_all <- plot_single_feature(
    feature = "Glucose",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_true("metab-t-clinical" %in% via_all$data$assay)

  # requesting another ome must not return clinical chemistry alongside it, and the
  # error has to name the ome that would have worked
  expect_error(plot_single_feature(feature = "Glucose",
                                   selected_omes = "transcript-rna-seq",
                                   selected_tissues = "blood",
                                   verbose = FALSE),
               "measured by clinical chemistry \\(metab-t-clinical\\)")

  # "metab" is the research platforms, and does not imply the clinical one. In blood the
  # analyte is also on the conventional panel, so the request is answerable without it —
  # what matters is that the clinical assay is not what comes back.
  research <- plot_single_feature(
    feature = "Glucose",
    selected_omes = "metab",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_false("metab-t-clinical" %in% research$data$assay)
  expect_true("metab-t-conv" %in% research$data$assay)

  # the clinical omes are separate from each other too
  expect_error(plot_single_feature(feature = "Glucose",
                                   selected_omes = "prot-clinical",
                                   selected_tissues = "blood",
                                   verbose = FALSE),
               "measured by clinical chemistry \\(metab-t-clinical\\)")

  # a feature that genuinely is not in the data keeps the original error
  expect_error(plot_single_feature(feature = "VEGFA",
                                   selected_omes = "metab-t-clinical",
                                   selected_tissues = "blood",
                                   verbose = FALSE),
               "No differential analysis corresponds")
})

test_that("the conventional metabolomics platform is plotted and labelled", {
  # metab-t-conv is no longer filtered out. Upstream assay_codes labels it "Conv(T)",
  # which does not distinguish it from its log2 twin metab-t-clinical, so
  # plot_single_feature() overrides the label rather than taking upstream verbatim.
  res <- plot_single_feature(
    feature = "Glucose",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_true("metab-t-conv" %in% res$data$assay)
  expect_false(any(is.na(res$data$tissue_assay)))
  expect_true("Blood Conv. Metab (log2)" %in% res$data$tissue_assay)
})

test_that("clinical prot analytes follow the same gate", {
  res <- plot_single_feature(
    feature = "Insulin",
    selected_omes = "prot-clinical",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_equal(unique(as.character(res$data$assay)), "prot-clinical")

  expect_error(plot_single_feature(feature = "Insulin",
                                   selected_omes = "prot-ol",
                                   selected_tissues = "blood",
                                   verbose = FALSE),
               "measured by clinical chemistry \\(prot-clinical\\)")
})

test_that("an analyte measured clinically and on a research platform splits by ome", {
  both <- plot_single_feature(
    feature = "Cortisol",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_true(all(c("metab-t-clinical", "metab-u-hilicpos") %in% both$data$assay))

  research_only <- plot_single_feature(
    feature = "Cortisol",
    selected_omes = "metab",
    selected_tissues = "blood",
    verbose = FALSE
  )
  expect_false("metab-t-clinical" %in% research_only$data$assay)
  expect_true("metab-u-hilicpos" %in% research_only$data$assay)
})

test_that("plot_single_feature rejects invalid tissue", {
  expect_error(
    plot_single_feature(
      feature = "VEGFA",
      selected_tissues = "liver"
    )
  )
})

test_that("plot_single_feature rejects invalid ome", {
  expect_error(
    plot_single_feature(
      feature = "VEGFA",
      selected_omes = "fake-ome"
    )
  )
})

test_that("plot_single_feature legends collect across tissues", {
  skip_if_not_installed("patchwork")

  tissues <- c("adipose", "blood", "muscle")
  plots <- lapply(tissues, function(tissue) {
    plot_single_feature(
      feature = "TAMALIN",
      selected_tissues = tissue,
      selected_omes = "transcript-rna-seq",
      verbose = FALSE
    )
  })

  # the reported failure needs a mix of tissues: at least one with no timepoint below
  # the p threshold, and at least one with some
  has_significant <- vapply(
    plots,
    function(p) any(p$data$below_p_cutoff == "Below p threshold"),
    logical(1)
  )
  skip_if_not(any(has_significant) && !all(has_significant),
              "TAMALIN no longer mixes significant and non-significant tissues")

  # every plot has to carry the same key set, or patchwork keeps one legend per variant
  fill_keys <- lapply(plots, ggplot2::get_guide_data, aesthetic = "fill")
  expect_equal(fill_keys[[2]], fill_keys[[1]])
  expect_equal(fill_keys[[3]], fill_keys[[1]])

  combined <- patchwork::wrap_plots(plots, ncol = 3, guides = "collect") &
    ggplot2::theme(legend.position = "right")
  assembled <- patchwork::patchworkGrob(combined)
  guide_box <- assembled$grobs[[which(assembled$layout$name == "guide-box")[1]]]

  # one legend for the exercise group colors, one for the p threshold fills
  expect_equal(sum(guide_box$layout$name == "guides"), 2)
})

test_that("plot_single_feature saves to file when output_file provided", {
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)

  res <- plot_single_feature(
    feature = "VEGFA",
    selected_tissues = "muscle",
    selected_omes = "transcript-rna-seq",
    output_file = tmp,
    verbose = FALSE
  )
  expect_true(file.exists(tmp))
  expect_s3_class(res, "ggplot")
})

# qc_data shaped like load_qc() output: vialLabel numeric, as it is read in, and Timepoint
# a factor on the study's levels
fake_qc_entry <- function(qc_norm, sample_metadata) {
  return(list(qc_norm = qc_norm, sample_metadata = sample_metadata,
              feature_metadata = data.frame()))
}
timepoint_levels <- c("pre_exercise", "during_20_min", "during_40_min", "post_10_min",
                      "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")

test_that(".summary_stats_from_qc summarises acute samples of the requested features", {
  qc_norm <- data.frame(`101` = c(1, 10), `102` = c(3, 20), `103` = c(5, 30), `104` = c(100, 100),
                        row.names = c("chr1:1-100", "chr1:200-300"), check.names = FALSE)
  sample_metadata <- data.frame(
    vialLabel = c(101, 102, 103, 104),
    visitcode = c("ADU_BAS", "ADU_BAS", "ADU_BAS", "ADU_POST"),
    randomGroupCode = c("ADUEndur", "ADUEndur", "ADUResist", "ADUEndur"),
    Timepoint = factor("pre_exercise", levels = timepoint_levels)
  )
  qc_data <- list(muscle = list(`epigen-atac-seq` = fake_qc_entry(qc_norm, sample_metadata)))
  features <- data.frame(tissue = "muscle", assay = "epigen-atac-seq",
                         feature_id = "chr1:1-100")

  res <- MotrpacHumanPreSuspensionAnalysis:::.summary_stats_from_qc(qc_data, features)

  expect_setequal(res$feature_id, "chr1:1-100")
  expect_named(res, c("tissue", "assay", "randomGroupCode", "Timepoint", "feature_id",
                      "Count", "Mean", "SD"), ignore.order = TRUE)
  endur <- res[res$randomGroupCode == "ADUEndur", ]
  # 104 is not an acute sample and is left out
  expect_equal(endur$Count, 2)
  expect_equal(endur$Mean, 2)
  expect_equal(endur$SD, sd(c(1, 3)))
  # one sample: the SD is NA, as it is in the shipped objects
  resist <- res[res$randomGroupCode == "ADUResist", ]
  expect_equal(resist$Count, 1)
  expect_true(is.na(resist$SD))

  # a tissue or ome qc_data does not carry, or a feature not in qc_norm, gives zero rows
  absent_ome <- data.frame(tissue = "blood", assay = "epigen-atac-seq", feature_id = "chr1:1-100")
  absent_feature <- data.frame(tissue = "muscle", assay = "epigen-atac-seq", feature_id = "chr9:1-2")
  expect_equal(nrow(MotrpacHumanPreSuspensionAnalysis:::.summary_stats_from_qc(qc_data, absent_ome)), 0)
  expect_equal(nrow(MotrpacHumanPreSuspensionAnalysis:::.summary_stats_from_qc(qc_data, absent_feature)), 0)
})

test_that(".summary_stats_from_qc handles several features, omes and tissues at once", {
  make_entry <- function(values, feature) {
    return(fake_qc_entry(
      data.frame(`201` = values[1], `202` = values[2], row.names = feature, check.names = FALSE),
      data.frame(vialLabel = c(201, 202), visitcode = "ADU_BAS", randomGroupCode = "ADUControl",
                 Timepoint = factor(c("pre_exercise", "post_24_hr"), levels = timepoint_levels))))
  }
  # tissue names are matched case-insensitively; metabolomics is keyed by platform
  qc_data <- list(MUSCLE = list(`epigen-atac-seq` = make_entry(c(1, 2), "chr1:1-100")),
                  blood = list(`metab-u-rppos` = make_entry(c(5, 7), "CAR 10:0")))
  features <- data.frame(tissue = c("muscle", "blood"),
                         assay = c("epigen-atac-seq", "metab-u-rppos"),
                         feature_id = c("chr1:1-100", "CAR 10:0"))

  res <- MotrpacHumanPreSuspensionAnalysis:::.summary_stats_from_qc(qc_data, features)

  expect_equal(nrow(res), 4)
  expect_setequal(paste(res$tissue, res$assay), c("muscle epigen-atac-seq", "blood metab-u-rppos"))
  expect_equal(res$Mean[res$assay == "metab-u-rppos" & res$Timepoint == "post_24_hr"], 7)
})

test_that(".summary_stats_from_qc drops duplicate, ungrouped and off-grid samples", {
  qc_norm <- data.frame(`301` = 1, `302` = 3, `303` = 50, `304` = 70,
                        row.names = "chr1:1-100", check.names = FALSE)
  sample_metadata <- data.frame(
    vialLabel = c(301, 302, 302, 303, 304),
    visitcode = "ADU_BAS",
    randomGroupCode = c("ADUEndur", "ADUEndur", "ADUEndur", NA, "ADUEndur"),
    Timepoint = c("pre_exercise", "pre_exercise", "pre_exercise", "pre_exercise", "not_a_timepoint")
  )
  qc_data <- list(muscle = list(`epigen-atac-seq` = fake_qc_entry(qc_norm, sample_metadata)))
  features <- data.frame(tissue = "muscle", assay = "epigen-atac-seq", feature_id = "chr1:1-100")

  expect_message(
    res <- MotrpacHumanPreSuspensionAnalysis:::.summary_stats_from_qc(
      qc_data, features, timepoints = timepoint_levels),
    "2 acute muscle epigen-atac-seq sample\\(s\\)")

  # 302 once, 303 without a group and 304 off the timepoint grid left out
  expect_equal(nrow(res), 1)
  expect_equal(res$Count, 2)
  expect_equal(res$Mean, 2)
})

test_that(".summary_stats_from_qc rejects qc_data that is not load_qc() output", {
  entry <- fake_qc_entry(data.frame(`401` = 1, row.names = "chr1:1-100", check.names = FALSE),
                         data.frame(vialLabel = 401, visitcode = "ADU_BAS",
                                    randomGroupCode = "ADUEndur", Timepoint = "pre_exercise"))
  features <- data.frame(tissue = "muscle", assay = "epigen-atac-seq", feature_id = "chr1:1-100")
  from_qc <- function(q) {
    return(MotrpacHumanPreSuspensionAnalysis:::.summary_stats_from_qc(q, features))
  }

  # one tissue, or one ome, of the load_qc() result instead of the whole of it
  expect_error(from_qc(list(`epigen-atac-seq` = entry)), "whole nested list")
  expect_error(from_qc(entry), "whole nested list")
  expect_error(from_qc(data.frame(x = 1)), "whole nested list")

  no_visit <- entry
  no_visit$sample_metadata$visitcode <- NULL
  expect_error(from_qc(list(muscle = list(`epigen-atac-seq` = no_visit))), "no visitcode column")

  not_numeric <- entry
  not_numeric$qc_norm[[1]] <- "1"
  expect_error(from_qc(list(muscle = list(`epigen-atac-seq` = not_numeric))), "not numeric")
})

# The epigenomic case in miniature, without the download: VEGFA is given differential
# analysis but no summary statistics, the way a non-significant ATAC peak is.
vegfa_id <- "ENSG00000112715.26"
local_sum_stats_without_vegfa <- function(env = parent.frame()) {
  real <- MotrpacHumanPreSuspensionAnalysis::load_summary_stats
  testthat::local_mocked_bindings(
    load_summary_stats = function(...) {
      out <- real(...)
      return(out[out$feature_id != vegfa_id, ])
    },
    .env = env)
  return(invisible(NULL))
}
vegfa_qc_data <- function() {
  shipped <- as.data.frame(MotrpacHumanPreSuspensionAnalysis::MUSCLE_TRANSCRIPT_RNA_SEQ_SUM_STATS)
  cells <- unique(shipped[shipped$feature_id == vegfa_id, c("randomGroupCode", "Timepoint")])
  # two vials per group and timepoint, valued so each mean is known: 10 * cell index
  vials <- seq_len(2 * nrow(cells)) + 500
  cell_index <- rep(seq_len(nrow(cells)), each = 2)
  qc_norm <- as.data.frame(t(setNames(10 * cell_index + c(-1, 1), vials)), check.names = FALSE)
  rownames(qc_norm) <- vegfa_id
  sample_metadata <- data.frame(vialLabel = vials, visitcode = "ADU_BAS",
                                randomGroupCode = cells$randomGroupCode[cell_index],
                                Timepoint = cells$Timepoint[cell_index])
  return(list(qc_data = list(muscle = list(`transcript-rna-seq` = fake_qc_entry(qc_norm, sample_metadata))),
              cells = cells))
}

test_that("plot_single_feature says when a feature has no summary statistics", {
  local_sum_stats_without_vegfa()
  expect_message(
    res <- plot_single_feature(feature = "VEGFA", selected_tissues = "muscle",
                               selected_omes = "transcript-rna-seq", verbose = FALSE),
    "No summary statistics for muscle transcript-rna-seq ENSG00000112715.26.*pass `qc_data`")
  expect_s3_class(res, "ggplot")
  expect_true(all(is.na(res$data$Mean)))
})

test_that("plot_single_feature computes missing summary statistics from qc_data", {
  local_sum_stats_without_vegfa()
  fake <- vegfa_qc_data()

  expect_no_message(
    res <- plot_single_feature(feature = "VEGFA", selected_tissues = "muscle",
                               selected_omes = "transcript-rna-seq", verbose = FALSE,
                               qc_data = fake$qc_data),
    message = "No summary statistics")

  expect_false(anyNA(res$data$Mean))
  expect_equal(sort(unique(res$data$Mean)), 10 * seq_len(nrow(fake$cells)))
  expect_true(all(res$data$Count == 2))
  # the x-axis keeps the study's timepoint order rather than sorting alphabetically
  expect_s3_class(res$data$Timepoint, "factor")
  plotted <- levels(droplevels(res$data$Timepoint))
  expect_equal(plotted[1], "Pre")
  expect_equal(plotted, intersect(c("Pre", "D20M", "D40M", "P10M", "P15-45M", "P3.5/4H", "P24H"),
                                  plotted))
})

test_that("plot_single_feature names what qc_data could not fill", {
  local_sum_stats_without_vegfa()
  fake <- vegfa_qc_data()
  rownames(fake$qc_data$muscle$`transcript-rna-seq`$qc_norm) <- "ENSG_SOMETHING_ELSE"

  expect_message(
    plot_single_feature(feature = "VEGFA", selected_tissues = "muscle",
                        selected_omes = "transcript-rna-seq", verbose = FALSE,
                        qc_data = fake$qc_data),
    "The `qc_data` you supplied does not carry these features")
})

test_that("qc_data leaves features with shipped summary statistics alone", {
  fake <- vegfa_qc_data()
  without <- plot_single_feature(feature = "VEGFA", selected_tissues = "muscle",
                                 selected_omes = "transcript-rna-seq", verbose = FALSE)
  with <- plot_single_feature(feature = "VEGFA", selected_tissues = "muscle",
                              selected_omes = "transcript-rna-seq", verbose = FALSE,
                              qc_data = fake$qc_data)
  expect_identical(with$data, without$data)
})

# Records every epigenomic load plot_single_feature asks for and stops there, so no
# file is downloaded; loads of the omes that ship in the package go through.
local_record_epigen_loads <- function(env = parent.frame()) {
  calls <- new.env()
  calls$epigen <- list()
  real <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis
  testthat::local_mocked_bindings(
    load_differential_analysis = function(..., epigen = FALSE) {
      args <- list(...)
      if (epigen) {
        calls$epigen[[length(calls$epigen) + 1]] <- args[c("selected_omes", "selected_tissues")]
        stop("epigenomic load recorded")
      }
      return(real(..., epigen = FALSE))
    },
    .env = env)
  return(calls)
}

test_that("epigen = TRUE downloads only the requested epigenomic tissue and ome", {
  calls <- local_record_epigen_loads()
  expect_error(plot_single_feature(feature = "chr15:84816704-84819087", epigen = TRUE,
                                   selected_tissues = "muscle",
                                   selected_omes = "epigen-atac-seq", verbose = FALSE),
               "epigenomic load recorded")
  expect_length(calls$epigen, 1)
  expect_equal(calls$epigen[[1]]$selected_omes, "epigen-atac-seq")
  expect_equal(calls$epigen[[1]]$selected_tissues, "muscle")
})

test_that("epigen = TRUE skips tissue and ome pairs that were not measured", {
  # adipose has methylcap but no ATAC
  calls <- local_record_epigen_loads()
  expect_error(plot_single_feature(feature = "ALPK3", epigen = TRUE,
                                   selected_tissues = c("adipose", "muscle"),
                                   selected_omes = "epigen-atac-seq", verbose = FALSE),
               "epigenomic load recorded")
  expect_equal(calls$epigen[[1]]$selected_tissues, "muscle")
})

test_that("epigen = TRUE downloads nothing when no epigenomic ome is requested", {
  calls <- local_record_epigen_loads()
  res <- plot_single_feature(feature = "VEGFA", epigen = TRUE, selected_tissues = "muscle",
                             selected_omes = "transcript-rna-seq", verbose = FALSE)
  expect_s3_class(res, "ggplot")
  expect_length(calls$epigen, 0)

  # adipose ATAC was never measured, so there is no file to fetch
  expect_error(plot_single_feature(feature = "ALPK3", epigen = TRUE, selected_tissues = "adipose",
                                   selected_omes = "epigen-atac-seq", verbose = FALSE),
               "No differential analysis corresponds")
  expect_length(calls$epigen, 0)
})

test_that("epigenomic omes requested with epigen = FALSE are reported as skipped", {
  expect_message(plot_single_feature(feature = "VEGFA", selected_tissues = "muscle",
                                     selected_omes = c("transcript-rna-seq", "epigen-atac-seq")),
                 "so epigenetic data will be skipped")
  expect_no_message(plot_single_feature(feature = "VEGFA", selected_tissues = "muscle",
                                        selected_omes = "transcript-rna-seq"),
                    message = "epigenetic data will be skipped")
})
