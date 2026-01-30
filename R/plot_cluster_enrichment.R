#' @title Bubble heatmap of CAMERA-PR or ORA results from fuzzy c-means
#'   clustering
#'
#' @description Generate a bubble heatmap from the results of nonparametric
#'   CAMERA-PR or ORA applied to the output of fuzzy c-means clustering.
#'
#' @param x a \code{data.frame}. Either
#'   \code{\link[MotrpacHumanPreSuspensionAnalysis]{FCM_CAMERA}} or
#'   \code{\link[MotrpacHumanPreSuspensionAnalysis]{FCM_ORA}}.
#' @param set_ids character or \code{NULL}; one or more of \code{x$set_id}. The
#'   dataset \code{x} will be filtered to these sets. If \code{NULL}, the
#'   \code{n_top} most significant sets from each contrast will be selected.
#' @param selected_tissues character; one or more tissues that will be used to
#'   create the heatmap. If multiple tissues are selected, only those terms
#'   tested in at least two tissues will be plotted.
#' @param selected_ome character; the ome that will be used to create the
#'   heatmap.
#' @param n_top integer; if \code{set_ids} is \code{NULL}, the \code{n_top} most
#'   significant molecular signatures from each tissue and cluster will be
#'   selected.
#' @param column_title character or \code{NULL}; title for the heatmap.
#' @param filename character; path to a PDF file to save the heatmap.
#' @param padj_cutoff numeric; cutoff in terms of p-value for an asterisk to be included
#'
#' @returns Nothing. The heatmap is saved to a PDF specified by \code{filename}.
#'
#' @import ComplexHeatmap
#' @importFrom dplyr %>% filter slice_min pull mutate arrange
#' @importFrom grDevices dev.off cairo_pdf
#' @importFrom grid convertUnit unit gpar
#' @importFrom latex2exp TeX
#'
#' @export plot_cluster_enrichment
#'
#' @author Tyler Sagendorf

plot_cluster_enrichment <- function(x,
                                    set_ids = NULL,
                                    selected_tissues = c("adipose",
                                                         "blood",
                                                         "muscle"),
                                    selected_ome = c("transcript-rna-seq",
                                                     "prot-pr",
                                                     "prot-ph",
                                                     "metab"),
                                    padj_cutoff = 0.05,
                                    n_top = 6L,
                                    column_title = NULL,
                                    filename)
{
  # Allow users with older R versions to still use the package
  check_package_installation(pkg = "TMSig", fun = "plot_cluster_enrichment")

  selected_tissues <- match.arg(selected_tissues,
                                choices = c("adipose", "blood", "muscle"),
                                several.ok = TRUE)

  selected_ome <- match.arg(selected_ome,
                            choices = c("transcript-rna-seq",
                                        "prot-pr",
                                        "prot-ph",
                                        "metab"))

  if (!is.character(filename) || !grepl("\\.pdf$", filename))
    stop("`filename` must be a path to a PDF file used to save the heatmap.")

  required_cols <- c("tissue", "assay", "cluster", "set_id", "set",
                     "set_short", "p_value", "adj_p_value")

  missing_cols <- setdiff(required_cols, colnames(x))

  if (length(missing_cols))
    stop("`x` is missing the following required column(s): ",
         paste(missing_cols, collapse = ", "))

  x <- x %>%
    filter(tissue %in% selected_tissues,
           assay == selected_ome) %>%
    droplevels.data.frame()

  if (nrow(x) == 0L)
    stop("Invalid combination of `selected_tissues` and `selected_ome`.")

  x <- x %>%
    filter(.by = set,
           any(adj_p_value < padj_cutoff))

  if (nrow(x) == 0L)
    stop("No sets significantly localize to any cluster at padj_cutoff = ",
         padj_cutoff)

  selected_tissues <- intersect(levels(x$tissue), selected_tissues)

  # Keep sets that appear in at least 2 tissues if multiple tissues are selected
  if (length(selected_tissues) > 1L) {
    x <- x %>%
      filter(.by = set,
             length(unique(tissue)) > 1L) %>%
      arrange(tissue, cluster) %>%
      mutate(cluster = paste(tissue, " ", cluster),
             cluster = factor(cluster, levels = unique(cluster)))
  }

  if (is.null(set_ids)) {
    # If set_ids is NULL, select the most significant sets from each contrast
    if (!is.numeric(n_top) && !is.infinite(n_top) && length(n_top) != 1L)
      stop("`n_top` must be a single integer or Inf.")

    n_top <- max(1L, n_top, na.rm = FALSE)

    set_ids <- x %>%
      filter(adj_p_value < padj_cutoff) %>%
      # p_value used for ordering to avoid ties
      slice_min(p_value,
                by = cluster,
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

  column_labels <- levels(x$cluster)

  if (length(unique(x$tissue)) > 1L) {
    column_labels <- sub(".* ", "", levels(x$cluster))

    # Separate columns by tissues
    column_split <- sub(" .*", "", levels(x$cluster)) %>%
      factor(levels = unique(.))

    # Top annotation for tissues
    anno_top <- HeatmapAnnotation(
      Tissue = column_split,
      col = list(Tissue = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[levels(column_split)]),
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
  } else {
    column_labels <- levels(x$cluster)
    column_split <- anno_top <- NULL
  }

  # Dynamic height and width
  height <- convertUnit(length(unique(x$set_id)) * unit(14, "pt"), "in") +
    unit(0.5, "in")
  height <- max(as.numeric(height), 4) + 0.3 + 0.25 * (!is.null(anno_top)) +
    0.25 * (!is.null(column_title))

  width <- convertUnit(nlevels(x$cluster) * unit(14, "pt"), "in") +
    convertUnit(max_text_width(unique(x$set_short),
                               gp = gpar(fontsize = 0.9 * 14)), "in")
  width <- as.numeric(width) + 2.5

  draw_args <- list()

  # Add padding below shorter heatmaps
  if (length(unique(x$set_id)) <= 10L) {
    draw_args <- list(
      padding = unit(c(100, 0, 0, 0), "pt")
    )
  }

  # Use cairo_pdf() to properly render ">=" symbol on Mac
  cairo_pdf(filename = filename,
            height = height,
            width = width,
            onefile = FALSE)

  x %>%
    mutate(log_p = -log10(p_value)) %>%
    TMSig::enrichmap(
      n_top = Inf,
      set_column = "set_short",
      statistic_column = "log_p",
      contrast_column = "cluster",
      padj_column = "adj_p_value",
      padj_legend_title = "BH Adjusted\nP-Value",
      colors = c("white", "#543483"),
      padj_cutoff = padj_cutoff,
      heatmap_args = list(
        na_col = "grey95",
        rect_gp = gpar(col = "grey85",
                       fill = "white"),
        column_names_side = "top",
        column_labels = column_labels,
        column_title_gp = gpar(fontsize = 12,
                               fontface = "bold"),
        column_title = column_title,
        column_split = column_split,
        top_annotation = anno_top,
        heatmap_legend_param = list(
          title = latex2exp::TeX("$\\bf{$-log$_{10}($P-Value$)}$")
        )
      ),
      draw_args = draw_args
    )

  invisible(dev.off())
}
