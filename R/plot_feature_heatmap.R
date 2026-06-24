#' @title Create a heatmap of select features or a choice pathway.
#'
#' @description Create a heatmap for user-specified features for a given tissue
#'   and ome combination. This function replaces `plot_pathway_features`
#'
#' @param feature_ids character or \code{NULL}; vector of feature IDs that will
#'   appear in the heatmap. Must be a subset of
#'   \code{HUMAN_FEATURE_TO_ID[["feature_id"]]}.
#' @param set_id character of \code{NULL}; if not \code{NULL}, the feature IDs
#'   in the set will be used to filter rows for the heatmap, but those IDs will
#'   not be used for the heatmap row labels.
#' @param platforms character or \code{NULL}; the platform(s) used to filter
#'   metabolites. If \code{NULL} (default), metabolites in \code{feature_ids}
#'   will be selected from all available platforms. If metabolites appear in
#'   more than one platform, the platform will appear before the metabolite name
#'   in the row names of the heatmap.
#' @param selected_tissue character; the tissue that will be used to create the
#'   heatmap.
#' @param selected_ome character; the ome that will be used to create the
#'   heatmap.
#' @param filename character; optional file name used to save the heatmap. If
#'   provided, the heatmap will not be drawn.
#' @param post_min numeric; for experiments, timepoints could be either 15, 30, or
#'   45 minutes depending on the analysis, this argument allows users to specify
#'   which of those 3 values to use, by default the value is NULL and will provide
#'   a generic label of post 15/30/45 min
#' @param post_hr numeric; for experiments, timepoints could be either 3.5 or 4 hours
#'   depending on the analysis, this argument allows users to specify
#'   which of those 2 values to use, by default the value is NULL and will provide
#'   a generic label of post 3.5/4 hr
#' @param full_modality_names logical; if TRUE the modality values are set as Endurance
#'   Exercise and Resistance Exercise but if FALSE the modality values are set as EE
#'   and RE (by default the value is FALSE)
#' @param contrast_type character; the type of contrasts to plot. One of
#'   "exercise_with_controls" (default), "exercise_no_controls",
#'   "Endur_vs_Resist", "baseline", or "control_only".
#' @param column_title character; the title you'd like to include for the columns. Usually empty
#' @param max_size numeric; largest number of pathways to display, if a pathway has too many features.
#'    Will automatically chose the first n pathways.
#' @param multi_tissue_clust_rows logical; if you'd like to cluster rows specifically in a multi-tissue case
#' @param verbose logical; for specific warnings and additional information.
#' @param ... Additional parameters to be added to a ComplexHeatmap call
#'
#' @returns Nothing. A heatmap is drawn or saved to a file if \code{filename} is
#'   provided.
#'
#' @export plot_feature_heatmap
#'
#' @author Tyler Sagendorf, Damon Leach, Christopher Jin
#'
#' @import ComplexHeatmap
#' @importFrom dplyr %>% filter mutate bind_rows
#' @importFrom grDevices dev.off
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider
#' @importFrom grid grid.text grid.rect convertUnit convertWidth convertHeight grobWidth grobHeight textGrob unit gpar
#' @examples
#' \dontrun{
#' plot_feature_heatmap(set_id = "11725",
#'   DA_list = DA_list,
#'   contrast_type = "exercise_with_controls",
#'   selected_tissue = "muscle",
#'   selected_ome = "prot-ph",
#'   filename = "sandbox/test_feature_heatmap.pdf")
#'}


plot_feature_heatmap <- function(feature_ids = NULL,
                                 set_id = NULL, # if provided, feature_ids is ignored
                                 platforms = NULL, # metabolomics only
                                 selected_tissue = c("adipose", "blood", "muscle"),
                                 selected_ome = c("transcript-rna-seq", "prot-pr",
                                                  "prot-ph", "metab", "prot-ol"),
                                 contrast_type = "exercise_with_controls",
                                 column_title = "",
                                 filename,
                                 max_size = NULL,
                                 post_min = NULL,
                                 post_hr = NULL,
                                 full_modality_names = FALSE,
                                 multi_tissue_clust_rows = FALSE,
                                 verbose = TRUE,
                                 ...)
{
  .validate_feature_heatmap_inputs(post_min,
                                   post_hr,
                                   selected_tissue,
                                   selected_ome,
                                   full_modality_names,
                                   contrast_type)

  DA_res = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(selected_omes = selected_ome,
                                      selected_tissues = selected_tissue,
                                      combine_with_featgene = TRUE,
                                      single_matrix = TRUE) %>%
    dplyr::filter(contrast_type == !!contrast_type) %>%
    dplyr::mutate(across(.cols = any_of(c("feature_id", "gene_symbol",
                                   "platform", "flanking_sequence")),
                  .fns = as.character))

  # If set_id is provided, select the genes/features from the set
  if (!is.null(set_id)) {
    if (!is.character(set_id) || length(set_id) != 1L) {
      stop("`set_id` must be a character string specifying a set ID ",
           "from the `SET_TO_ID` object. Make sure to include preceding 0s if there are any")
    }
    pathway <- MotrpacHumanPreSuspensionAnalysis::SET_TO_ID %>%
      dplyr::filter(set_id == !!set_id) %>%
      dplyr::pull(set)

    if(length(pathway) == 0) stop("`set_id` must be a character string specifying a set ID ",
                                  "from the `SET_TO_ID` object. Make sure to include preceding 0s if there are any")

    pathways <- unlist(structure(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES, names = NULL),
                       recursive = FALSE)
    feature_ids <- pathways[[pathway]]
    #subset to first n features
    if(!is.null(max_size))
      feature_ids = feature_ids[1:max_size]
  }

  if (!is.vector(feature_ids, mode = "character")) {
    stop("If `set_id` is not provided, `feature_ids` must be a character ",
         "vector of feature IDs that will appear in the heatmap.")
  }

  if (is.null(set_id)) {
    # Filter before modifying IDs
    DA_res <- dplyr::filter(DA_res, feature_id %in% feature_ids)
  }

  if (!nrow(DA_res)) {
    stop("`feature_ids` do not match feature IDs in the ",
         sprintf("%s %s differential analysis results.",
                 selected_tissue, selected_ome))
  }

  if (selected_ome %in% c("prot-pr", "transcript-rna-seq", "prot-ol")) {
    DA_res <- DA_res %>%
      dplyr::mutate(feature_id = ifelse(!is.na(gene_symbol),
                                 gene_symbol,
                                 feature_id))
  } else if (selected_ome == "prot-ph") {
    DA_res <- DA_res %>%
      dplyr::mutate(feature_id = gsub("[sty]", ";", feature_id),
             feature_id = sub("(.*);$", "\\1", feature_id),
             feature_id = sub(".*_", "", feature_id),
             feature_id = paste(gene_symbol, feature_id))

  } else  if (selected_ome == "metab" && !is.null(platforms)) {
    platforms <- match.arg(arg = platforms,
                           choices = paste0("metab-",
                                            c("t-amines", "t-conv", "t-imm-crt",
                                              "t-imm-glc", "t-imm-ins",
                                              "t-oxylipneg", "t-tca",
                                              "u-hilicpos", "u-ionpneg",
                                              "u-lrpneg", "u-lrppos",
                                              "u-rpneg", "u-rppos")),
                           several.ok = TRUE)

    DA_res <- DA_res %>%
      dplyr::filter(platform %in% platforms)

    if (!nrow(DA_res)) {
      stop("Metabolites in `feature_ids` not found in the provided platforms.")
    }
  }

  if (!is.null(set_id)) {
    if (selected_ome %in% c("prot-pr", "transcript-rna-seq")) {
      DA_res <- dplyr::filter(DA_res, feature_id %in% feature_ids)
    } else if (selected_ome == "prot-ph") {
      DA_res <- dplyr::filter(DA_res, gene_symbol %in% feature_ids)
    } else if (selected_ome == "metab") {
      DA_res <- dplyr::filter(DA_res, feature_id %in% feature_ids)
    }
  }

  if (nrow(DA_res) == 0L) {
    stop("The combination of `set_id`, `selected_tissue`, and `selected_ome` ",
         "is not valid.")
  }


  x <- DA_res %>%
    dplyr::mutate(contrast2 = paste(tissue, contrast))

  # get factor level in order for contrast2
  contrast2_levels_df = x %>%
    dplyr::mutate(Timepoint = factor(Timepoint, levels = c("during_20_min","during_40_min","post_10_min",
                                                           "post_15_30_45_min","post_3.5_4_hr","post_24_hr"))) %>%
    dplyr::arrange(tissue,randomGroupCode,Timepoint) %>%
    dplyr::select(tissue,randomGroupCode,Timepoint,contrast2) %>%
    dplyr::distinct()
  contrast2_levels <- contrast2_levels_df$contrast2

  # updated x
  x <- x %>%
    dplyr::mutate(contrast2 = factor(contrast2, levels = contrast2_levels)) %>%
    dplyr::arrange(contrast2, abs(z.std), feature_id) %>%
    dplyr::filter(.by = contrast2,
           !duplicated(feature_id)) %>%
    dplyr::select(feature_id, contrast, contrast2, z.std, adj_p_value,tissue) %>%
    droplevels.data.frame()

  # Better contrast labels
  contrast_df <- .add_contrast_labels() %>%
    dplyr::filter(contrast %in% levels(x$contrast)) %>%
    droplevels.data.frame()

  contrast_colors_vector <- .contrast_colors()

  if(!is.null(post_min)){
    post_min_str = paste0("post ", post_min, " min")
    # update level information
    new_levels = levels(contrast_df$contrast_labels)
    new_levels[new_levels == "post 15/30/45 min"] = post_min_str
    # update values in df
    contrast_df$contrast_labels <- as.character(contrast_df$contrast_labels)
    contrast_df$contrast_labels[contrast_df$contrast_labels == "post 15/30/45 min"] = post_min_str
    # reset to factor with new levels
    contrast_df$contrast_labels <- factor(contrast_df$contrast_labels, levels = new_levels)
    # update color vectors
    names(contrast_colors_vector)[names(contrast_colors_vector) == "post 15/30/45 min"] = post_min_str
  }

  if(!is.null(post_hr)){
    post_hr_str = paste0("post ", post_hr, " hr")
    # update level information
    new_levels = levels(contrast_df$contrast_labels)
    new_levels[new_levels == "post 3.5/4 hr"] = post_hr_str
    # update values in df
    contrast_df$contrast_labels <- as.character(contrast_df$contrast_labels)
    contrast_df$contrast_labels[contrast_df$contrast_labels == "post 3.5/4 hr"] = post_hr_str
    # reset to factor with new levels
    contrast_df$contrast_labels <- factor(contrast_df$contrast_labels, levels = new_levels)
    # update color vectors
    names(contrast_colors_vector)[names(contrast_colors_vector) == "post 3.5/4 hr"] = post_hr_str
  }

  ## Create heatmap ----
  column_df <- distinct(x,tissue,contrast) %>%
    dplyr::mutate(tissue = factor(tissue, levels = selected_tissue)) %>%
    dplyr::arrange(contrast,tissue) %>%
    dplyr::left_join(contrast_df,
              by = "contrast") %>%
    dplyr::rename(modality = anno_group) %>%
    dplyr::mutate(modality = factor(modality,
                             levels = c("EE", "RE")))

  if(full_modality_names == TRUE){
    column_df$modality <- as.character(column_df$modality)
    column_df$modality[column_df$modality == "EE"] = "Endurance Exercise"
    column_df$modality[column_df$modality == "RE"] = "Resistance Exercise"
    column_df$modality <- factor(column_df$modality, levels = c("Endurance Exercise","Resistance Exercise"))
  }

  # order the annotations by tissue -> modality -> contrast_labels
  anno_df <- dplyr::select(column_df,
                           Tissue = tissue,
                           Modality = modality,
                           Timepoint = contrast_labels) %>%
    dplyr::arrange(Tissue,Modality,Timepoint)

  if(full_modality_names == TRUE){
    anno_col <- list(
      "Tissue" = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[selected_tissue],
      "Modality" = setNames(c("#d95f02", "#1b9e77"),
                            c("Endurance Exercise", "Resistance Exercise")),
      "Timepoint" = contrast_colors_vector[levels(anno_df$Timepoint)]
    )
  } else {
    anno_col <- list(
      "Tissue" = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[selected_tissue],
      "Modality" = setNames(c("#d95f02", "#1b9e77"),
                            c("EE", "RE")),
      "Timepoint" = contrast_colors_vector[levels(anno_df$Timepoint)]
    )
  }


  # contrasts from different tissues must be treated as distinct
  column_df <- column_df %>%
    dplyr::mutate(contrast2 = paste(tissue, contrast),
                  contrast2 = factor(contrast2, levels = unique(contrast2)))
  # If there is a single tissue, remove the tissue annotation
  if (length(selected_tissue) == 1L) {
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

  if (length(anno_col) > 0 ) {
    if (!is.null(anno_df[["Tissue"]]) && !is.null(anno_df[["Modality"]])) {
      column_split <- anno_df %>%
        dplyr::mutate(column_split = paste(Tissue, Modality),
               row_order = 1:n()) %>%
        dplyr::arrange(Tissue, Modality) %>%
        dplyr::mutate(column_split = factor(column_split,
                                     levels = unique(column_split))) %>%
        dplyr::arrange(row_order) %>%
        dplyr::pull(column_split)
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
  } else {
    top_annotation <- column_split <- NULL
  }

  n_features <- length(unique(x[["feature_id"]]))

  # Dynamic heatmap height. Each cell of the heatmap is 14 pt, by default
  height <- convertUnit(n_features * unit(14, "pt"), "in")
  height <- max(as.numeric(height), 4.5) + 1.5 + show_column_names

  # Dynamic heatmap width
  row_label_width <- ComplexHeatmap::max_text_width(
    text = x[["feature_id"]],
    gp = gpar(fontsize = 0.9 * 14)
  )
  row_label_width <- convertUnit(row_label_width, "in")

  width <- convertUnit(nlevels(x[["contrast"]]) * unit(14, "pt"), "in")
  width_extra = ifelse(length(selected_tissue) > 1, 5, 3)
  width <- as.numeric(width + row_label_width) + width_extra


  draw_args <- list(
    annotation_legend_list = NULL,
    merge_legends = TRUE
  )

  if (n_features <= 10L) {
    draw_args[["padding"]] <- unit(c(80, 0, 0, 0), "pt")
  }

  extended_range <- TMSig::extendRangeNum(x[["z.std"]], nearest = 0.1)

  if (all(extended_range <= 0)) {
    breaks <- c(extended_range[1], 0)
  } else if (all(extended_range >= 0)) {
    breaks <- c(0, extended_range[2])
  } else {
    max_r <- max(abs(extended_range))
    breaks <- c(extended_range[1], 0, extended_range[2])
  }

  # if we are comparing multiple tissues, likely to have
  # missing values which likely will cause issues with clustering
  # if we are comparing multiple tissues, likely to have
  # missing values which likely will cause issues with clustering
  # thus, either filter features with missing values, or don't cluster rows
  if (multi_tissue_clust_rows == TRUE){
    x<-x %>%
      dplyr::filter(feature_id %in% test_x$feature_id)
    clust_row_info == TRUE
  } else {
    clust_row_info <- ifelse(length(selected_tissue) > 1, FALSE, TRUE)
  }

  TMSig::enrichmap(
    x = x,
    n_top = Inf,
    set_column = "feature_id",
    statistic_column = "z.std",
    contrast_column = "contrast2",
    padj_column = "adj_p_value",
    plot_sig_only = FALSE,
    filename = filename,
    height = height,
    width = width,
    heatmap_color_fun = .feature_color_function,
    heatmap_args = list(
      layer_fun = .feature_layer_fun,
      cluster_rows = clust_row_info,
      column_split = column_split,
      column_labels = contrast_df$contrast_labels,
      show_column_names = show_column_names,
      column_names_side = "top",
      na_col = "grey80",
      column_title = gt_render(column_title,
                               padding = unit(c(0, 0, 0, 0.8), "in")),
      column_title_gp = gpar(fontsize = 12),
      top_annotation = top_annotation,
      heatmap_legend_param = list(
        title = "Z-Score",
        at = breaks,
        labels = breaks
      )
    ),

    draw_args = draw_args
  )
}


## Helper functions ------------------------------------------------------------

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

#function to make the names::color mapping to the aesthetic desired.
.contrast_colors = function() {
  tp_mod_colors =
    MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS %>%
    { setNames(., case_when(
      names(.) %in% c("during_20_min", "During", "during 20 min") ~ "during 20 min",
      names(.) %in% c("during_40_min", "during 40 min") ~ "during 40 min",
      names(.) %in% c("post_10_min", "post 10 min") ~ "post 10 min",
      names(.) %in% c("post_15_min", "post_30_min", "post_45_min", "post_15_30_45_min", "Early") ~ "post 15/30/45 min",
      names(.) %in% c("post_3.5_4_hr", "post_4_hr", "Mid") ~ "post 3.5/4 hr",
      names(.) %in% c("post_24_hr", "Late") ~ "post 24 hr",
      TRUE ~ names(.)  # keep any unmatched names unchanged
    )) }
  tp_mod_colors

}

.feature_layer_fun <- function(j, i, x, y, w, h, f) {
  grid.rect(x = x, y = y, width = w, height = h,
            gp = gpar(col = heatmap_args[["rect_gp"]][["col"]],
                      fill = f)
  )
  gb = textGrob("*")
  gb_w = convertWidth(grobWidth(gb), "mm")
  gb_h = convertHeight(grobHeight(gb), "mm")
  grid.text("*",
            x = x,
            y = y - gb_h*0.5 + gb_w*0.4,
            #r = pindex(dmat, i, j) / 2 * cell_size,
            # Significant bubbles get a black outline to separate from padj_fill
            gp = gpar(col = ifelse(pindex(padj_mat, i, j) < padj_cutoff,
                                   "black", NA))
  )
}


.feature_color_function <- function(statistics,
                                    colors = c("#3366ff", "darkred"))
{
  r <- range(statistics, na.rm = TRUE)

  # Extend range of values out to the nearest tenth
  extended_range <- TMSig::extendRangeNum(r, nearest = 0.1)

  max_r <- max(abs(extended_range))

  if (all(r >= 0)) {
    breaks <- c(0, +1) * max_r
    colors <- c("white", colors[2])
  } else if (all(r <= 0)) {
    breaks <- c(-1, 0) * max_r
    colors <- c(colors[1], "white")
  } else {
    breaks <- c(-1, 0, 1) * max_r
    colors <- c(colors[1], "white", colors[2])
  }

  return(list(breaks = breaks, colors = colors))
}

.validate_feature_heatmap_inputs = function(post_min,
                                            post_hr,
                                            selected_tissue,
                                            selected_ome,
                                            full_modality_names,
                                            contrast_type){
  if(!is.null(post_min)){
    if(length(post_min) > 1){
      stop("If not NULL, post_min can only have length 1")
    }
    if(!is.numeric(post_min)){
      stop("If not NULL, post_min must be numeric")
    }
    if(!post_min %in% c(15,30,45)){
      stop("If not NULL, post_min can only take on the numeric values of 15, 30, or 45")
    }
  }

  # check that post_hr is formatted correctly
  if(!is.null(post_hr)){
    if(length(post_hr) > 1){
      stop("If not NULL, post_hr can only have length 1")
    }
    if(!is.numeric(post_hr)){
      stop("If not NULL, post_hr must be numeric")
    }
    if(!post_hr %in% c(3.5,4)){
      stop("If not NULL, post_hr can only take on the numeric values of 3.5 or 4")
    }
  }

  # check that full_modality_names is formatted correctly
  if(length(full_modality_names) > 1){
    stop("full_modality_names can only have length 1")
  }
  if(!is.logical(full_modality_names)){
    stop("full_modality_names must be a logical argument")
  }


  selected_tissue <- match.arg(selected_tissue,
                               choices = c("adipose", "blood", "muscle"),
                               several.ok = TRUE)

  if(length(selected_ome) > 1)
    stop("only 1 ome at a time is supported. Please specify selected_ome input")

  selected_ome <- match.arg(selected_ome,
                            choices = c("transcript-rna-seq", "prot-pr",
                                        "prot-ph", "metab", "prot-ol"),
                            several.ok = FALSE)

  contrast_type <- match.arg(contrast_type,
                             choices = c("exercise_with_controls",
                                         "exercise_no_controls",
                                         "Endur_vs_Resist",
                                         "baseline",
                                         "control_only"),
                             several.ok = FALSE)
}

