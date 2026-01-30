#' @title Create a bubble heatmap of molecular signatures from CAMERA-PR or
#'   PTM-SEA results
#'
#' @description Create a heatmap of molecular signatures for a given tissue,
#'   ome, and contrast type combination. Only works with CAMERA-PR or PTM-SEA
#'   results.
#'
#' @param x a \code{data.frame} of CAMERA-PR or PTM-SEA results.
#' @param set_ids character or \code{NULL}; one or more of \code{x$set_id}. The
#'   dataset \code{x} will be filtered to these sets. If \code{NULL}, the top 6
#'   most significant sets from each contrast will be selected.
#' @param selected_tissues character; one or more tissues that will be used to
#'   create the heatmap. If multiple tissues are selected, only those terms
#'   tested in at least two tissues will be plotted.
#' @param selected_ome character; the ome that will be used to create the
#'   heatmap.
#' @param contrast_type character; the type of contrasts to plot. One of
#'   "exercise_with_controls" (default), "exercise_no_controls",
#'   "Endur_vs_Resist", "baseline", or "control_only".
#' @param padj_cutoff numeric; the p-value cutoff to add an asterisk indicating statistical significance.
#' @param n_top integer; if \code{set_ids} is \code{NULL}, the \code{n_top} most
#'   significant molecular signatures from each tissue and contrast will be
#'   selected.
#' @param zscore_colors character; length 2 vector of colors for the smallest
#'   and largest z-scores. Default is "#3366ff" (blue) and "darkred".
#' @param filename character; path to a PDF file to save the heatmap.
#'
#' @returns Nothing. The heatmap is saved to a PDF file.
#'
#' @author Tyler Sagendorf
#'
#' @export plot_enrich_heatmap
#'
#' @import ComplexHeatmap
#' @importFrom dplyr %>% filter slice_min pull distinct mutate left_join arrange
#'   rename count
#' @importFrom grDevices dev.off cairo_pdf
#' @importFrom grid convertUnit unit gpar
#' @importFrom latex2exp TeX

plot_enrich_heatmap <- function(x,
                                set_ids = NULL,
                                selected_tissues = c("adipose",
                                                     "blood",
                                                     "muscle"),
                                selected_ome = c("transcript-rna-seq",
                                                 "prot-pr",
                                                 "prot-ph",
                                                 "metab",
                                                 "prot-ol"),
                                contrast_type = c("exercise_with_controls",
                                                  "exercise_no_controls",
                                                  "Endur_vs_Resist",
                                                  "baseline",
                                                  "control_only"),
                                padj_cutoff = 0.05,
                                n_top = 6L,
                                zscore_colors = c("#3366ff", "darkred"),
                                filename)
{
  # Allow users with older R versions to still use the package
  check_package_installation(pkg = "TMSig", fun = "plot_bubble_heatmap")

  selected_tissues <- match.arg(selected_tissues,
                                choices = c("adipose", "blood", "muscle"),
                                several.ok = TRUE)

  selected_ome <- match.arg(selected_ome,
                            choices = c("transcript-rna-seq", "prot-pr",
                                        "prot-ph", "metab", "prot-ol"))

  contrast_type <- match.arg(contrast_type,
                             choices = c("exercise_with_controls",
                                         "exercise_no_controls",
                                         "Endur_vs_Resist",
                                         "baseline",
                                         "control_only"))

  ## Prepare CAMERA-PR results ----
  required_cols <- c("tissue", "assay", "contrast", "contrast_type",
                     "set_id", "set", "set_short", "p_value", "adj_p_value")

  missing_cols <- setdiff(required_cols, colnames(x))

  if (length(missing_cols))
    stop("`x` is missing the following required column(s): ",
         paste(missing_cols, collapse = ", "))

  if (!any(c("z.std", "NES") %in% colnames(x)))
    stop("`x` must include a 'z.std' or 'NES' column, depending on ",
         "whether it was produced by run_cameraPR or run_PTMSEA, ",
         "respectively.")

  x <- x %>%
    filter(contrast_type == !!contrast_type,
           assay == selected_ome,
           tissue %in% selected_tissues) %>%
    droplevels.data.frame()

  selected_tissues <- intersect(selected_tissues,
                                unique(x$tissue))

  # Keep sets that appear in at least 2 tissues if multiple tissues are selected
  if (length(selected_tissues) > 1L) {
    x <- x %>%
      filter(.by = set,
             length(unique(tissue)) > 1L)
  }

  x <- x %>%
    filter(.by = set,
           any(adj_p_value < padj_cutoff))

  if (!nrow(x))
    stop("No sets are significant at padj_cutoff = ", padj_cutoff)

  if (is.null(set_ids)) {
    # If set_ids is NULL, select the most significant sets from each contrast
    if (!is.numeric(n_top) && !is.infinite(n_top) && length(n_top) != 1L)
      stop("`n_top` must be a single integer or Inf.")

    n_top <- max(1L, n_top, na.rm = FALSE)

    set_ids <- x %>%
      filter(adj_p_value < padj_cutoff) %>%
      # p_value used for ordering to avoid ties
      slice_min(p_value,
                by = c(tissue, contrast),
                n = n_top) %>%
      pull(set_id) %>%
      unique() %>%
      as.character()
  } else {
    if (!is.character(set_ids))
      stop("`set_ids` must be NULL or a character vector of set identifiers ",
           "selected from SET_TO_ID$set_id")

    # Strip away any non digit characters (e.g. special formatting chars)
    set_ids <- gsub("[^[:digit:]]", "", set_ids)
    set_ids <- unique(set_ids)
    set_ids <- as.numeric(set_ids)
    set_ids <- sprintf("%05d", set_ids)

    if (all(!set_ids %in% x$set_id)) {
      stop("None of the `set_ids` are valid. Check `set_ids` and the values ",
           "of `selected_ome`, `selected_tissues`, and `padj_cutoff`.")
    }
  }

  x <- filter(x, set_id %in% set_ids)

  # Better contrast labels
  contrast_df <- .add_contrast_labels() %>%
    filter(contrast %in% levels(x$contrast)) %>%
    droplevels.data.frame()

  # If FALSE, assume PTM-SEA results
  is_camera <- !is.null(x[["z.std"]])

  ## Create heatmap ----
  column_df <- distinct(x, tissue, contrast) %>%
    mutate(tissue = factor(tissue, levels = selected_tissues)) %>%
    arrange(contrast, tissue) %>%
    left_join(contrast_df,
              by = "contrast") %>%
    rename(modality = anno_group) %>%
    mutate(modality = factor(modality,
                             levels = c("EE", "RE")))

  anno_df <- select(column_df,
                    Tissue = tissue,
                    Modality = modality,
                    Timepoint = contrast_labels)

  anno_col <- list(
    "Tissue" = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[selected_tissues],
    "Modality" = setNames(c("#d95f02", "#1b9e77"),
                          c("EE", "RE")),
    "Timepoint" = .contrast_colors()[levels(anno_df$Timepoint)]
  )

  # Contrasts from different tissues must be treated as distinct
  column_df <- column_df %>%
    mutate(contrast2 = paste(tissue, contrast),
           contrast2 = factor(contrast2, levels = unique(contrast2)))

  # If there is a single tissue, remove the tissue annotation
  if (length(selected_tissues) == 1L) {
    anno_df$Tissue <- NULL
    anno_col["Tissue"] <- NULL
  }

  show_column_names <- FALSE

  if (contrast_type == "baseline") {
    anno_df$Timepoint <- NULL
    anno_col["Timepoint"] <- NULL

    show_column_names <- TRUE
  }

  if (!contrast_type %in% c("exercise_with_controls", "exercise_no_controls")) {
    anno_df$Modality <- NULL
    anno_col["Modality"] <- NULL
  }

  top_annotation <- NULL
  column_split <- NULL

  if (length(anno_col)) {
    if (!is.null(anno_df[["Tissue"]]) && !is.null(anno_df[["Modality"]])) {
      column_split <- anno_df %>%
        mutate(column_split = paste(Tissue, Modality),
               row_order = 1:n()) %>%
        arrange(Tissue, Modality) %>%
        mutate(column_split = factor(column_split,
                                     levels = unique(column_split))) %>%
        arrange(row_order) %>%
        pull(column_split)
    } else if (!is.null(anno_df[["Tissue"]])) {
      column_split <- anno_df[["Tissue"]]
    } else if (!is.null(anno_df[["Modality"]])) {
      column_split <- anno_df[["Modality"]]
    }

    top_annotation <- HeatmapAnnotation(
      df = anno_df,
      col = anno_col,
      which = "column",
      border = TRUE,
      gap = unit(2, "pt"),
      annotation_name_gp = gpar(fontsize = 0.9 * 14),
      annotation_legend_param = list(
        border = TRUE,
        title_gp = gpar(fontsize = 0.9 * unit(14, "pt"),
                        fontface = "bold"),
        labels_gp = gpar(fontsize = 0.9 * unit(14, "pt"))
      )
    )
  }

  n_sets <- length(unique(x[["set"]]))

  # Dynamic heatmap height. Each cell of the heatmap is 14 pt, by default
  height <- convertUnit(n_sets * unit(14, "pt"), "in")
  height <- max(as.numeric(height), 4 + 2 * (1L - show_column_names)) +
    3 + show_column_names

  # Dynamic heatmap width
  row_label_width <- ComplexHeatmap::max_text_width(
    text = x[["set_short"]],
    gp = gpar(fontsize = 0.9 * 14)
  )
  row_label_width <- convertUnit(row_label_width, "in")

  # Dynamic heatmap width. Each cell of the heatmap is 14 pt, by default
  width <- convertUnit(nrow(anno_df) * unit(14, "pt"), "in")
  width <- as.numeric(width + row_label_width) + 3

  if (is_camera) {
    extended_range <- TMSig::extendRangeNum(x[["z.std"]], nearest = 0.1)
    heatmap_title <- "Z-Score"
    statistic_column <- "z.std"
  } else {
    extended_range <- TMSig::extendRangeNum(x[["NES"]], nearest = 0.1)
    heatmap_title <- statistic_column <- "NES"
  }

  # Breaks for heatmap color legend
  if (all(extended_range <= 0)) {
    breaks <- c(extended_range[1], 0)
  } else if (all(extended_range >= 0)) {
    breaks <- c(0, extended_range[2])
  } else {
    max_r <- max(abs(extended_range))
    breaks <- c(extended_range[1], 0, extended_range[2])
  }

  on.exit(invisible(dev.off())) # close device after drawing heatmap

  # Use cairo_pdf so that the ">=" symbol renders properly on Mac
  grDevices::cairo_pdf(
    filename = filename, height = height, width = width,
    onefile = FALSE, fallback_resolution = 300
  )

  # If at least one set is not measured in >= 50% of contrasts across tissues,
  # turn off row clustering because there is a chance it will fail anyway.
  cluster_rows <- x %>%
    count(set_short) %>% # number of nonmissing values
    pull(n) %>%
    {all(. >= 0.5 * max(.))}

  x <- x %>%
    mutate(contrast2 = paste(tissue, contrast),
           contrast2 = factor(contrast2, levels = levels(column_df$contrast2))) %>%
    mutate(set_short = factor(set_short, levels = sort(unique(set_short))))

  # If there are few terms, add padding to the bottom of the heatmap to avoid
  # increasing overall height of the file
  draw_args <- list()

  if (length(unique(x$set)) <= 10L) {
    draw_args <- list(
      padding = unit(c(185, 0, 0, 0), "pt")
    )
  } else if (length(unique(x$set)) <= 20L) {
    draw_args <- list(
      padding = unit(c(70, 0, 0, 0), "pt")
    )
  }

  # Bubble heatmap
  TMSig::enrichmap(
    x = x,
    n_top = Inf,
    set_column = "set_short",
    statistic_column = statistic_column,
    contrast_column = "contrast2",
    padj_column = "adj_p_value",
    padj_cutoff = padj_cutoff,
    plot_sig_only = TRUE,
    heatmap_color_fun = .enrich_heatmap_color_function,
    colors = zscore_colors,
    padj_legend_title = "BH Adjusted\nP-Value",
    draw_args = draw_args,
    heatmap_args = list(
      na_col = "grey95",
      rect_gp = gpar(fill = "white",
                     col = "grey85"),
      cluster_rows = cluster_rows,
      column_split = column_split,
      row_labels = latex2exp::TeX(levels(x$set_short)),
      column_labels = column_df$contrast_labels,
      show_column_names = show_column_names,
      column_names_side = "top",
      column_title_gp = gpar(fontsize = 0),
      top_annotation = top_annotation,
      heatmap_legend_param = list(
        title = heatmap_title,
        at = breaks,
        labels = breaks
      )
    )
  )
}



## Helper functions ------------------------------------------------------------

#' @title Create better contrast labels
#'
#' @returns A \code{data.frame} with all columns from \code{CONTRAST_CONVERTER}
#'   in addition to colums "contrast_label" (better contrast labels) and
#'   "anno_group" (either "Endur" or "Resist", only for exercise comparisons).
#'
#' @noRd

.add_contrast_labels <- function() {
  contrast_labels <- c(
    # exercise_with_controls and exercise_no_controls
    rep(c("during 20 min", "during 40 min", "post 10 min", "post 15/30/45 min",
          "post 3.5/4 hr", "post 24 hr", "post 10 min", "post 15/30/45 min",
          "post 3.5/4 hr", "post 24 hr"), times = 2L),
    # Endur vs. Resist
    c("post 10 min", "post 15/30/45 min", "post 3.5/4 hr", "post 24 hr"),
    # baseline
    c("Endur - Resist", "Endur - Control", "Resist - Control"),
    # control_only
    c("during 20 min", "during 40 min", "post 10 min", "post 15/30/45 min",
      "post 3.5/4 hr", "post 24 hr")
  )

  contrast_labels <- factor(
    x = contrast_labels,
    levels = c("Endur - Resist", "Endur - Control", "Resist - Control",
               "during 20 min", "during 40 min", "post 10 min",
               "post 15/30/45 min", "post 3.5/4 hr", "post 24 hr")
  )

  anno_group <- c(
    # exercise_with_controls and exercise_no_controls
    rep(rep(c("EE", "RE"), c(6L, 4L)), times = 2L),
    # Endur_vs_Resist and baseline, and control_only
    rep(NA, 13L)
  )

  out <- cbind(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER, contrast_labels, anno_group)

  return(out)
}


#' @title enrich_heatmap color function
#'
#' @param statistics numeric vector or matrix of enrichment analysis statistics
#'   (z-scores for CAMERA-PR or NES for PTM-SEA).
#'
#' @noRd

.enrich_heatmap_color_function <- function(statistics,
                                           colors = c("#3366ff", "darkred"))
{
  r <- range(statistics, na.rm = TRUE)

  # Extend range of values out to the nearest tenth
  extended_range <- TMSig::extendRangeNum(r, nearest = 0.1)

  if (all(r >= 0)) {
    breaks <- c(0, extended_range[2L])
    colors <- c("white", colors[2L])
  } else if (all(r <= 0)) {
    breaks <- c(extended_range[1L], 0)
    colors <- c(colors[1L], "white")
  } else {
    max_r <- max(abs(extended_range))
    breaks <- c(-max_r, 0, max_r)
    colors <- c(colors[1L], "white", colors[2L])
  }

  return(list(breaks = breaks, colors = colors))
}


#' @title Contrast colors
#'
#' @description Assign colors to select contrasts for heatmaps.
#'
#' @noRd
.contrast_colors <- function() {
  structure(
    c("#fde725",
      "#bad071",
      "#d1bbd7",
      "#ae76a3",
      "#882e72",
      "#61194f"),
    names = c("during 20 min",
              "during 40 min",
              "post 10 min",
              "post 15/30/45 min",
              "post 3.5/4 hr",
              "post 24 hr")
  )
}
