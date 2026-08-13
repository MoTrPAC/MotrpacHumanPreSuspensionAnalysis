#' @title Single feature plot function
#' @description Creates a single plot or faceted plot for a given feature or gene
#' across the selected tissues. Single plots are optimized for H = 1.54 in and W = 1.225 inches.
#' If a right hand legend is added, the width should be 1.715. This can be toggled
#' as needed by saving output and editing as a ggplot object
#'
#' @param feature A single character string corresponding to a feature_id, gene_symbol, or refmet_name within the human feature-to-gene table
#' @param p_level Numeric threshold for adjusted p value significance in differential analysis (so this would make individual points highlighted in black if below this threshold)
#' @param selected_tissues character; one of tissue_available_list.
#' @param selected_omes character; one of ome_available_list. The clinical
#'   chemistry omes in \code{\link{clinical_ome_list}()} are plotted only when
#'   they are requested, either by name or through \code{"all"}; asking for
#'   another ome no longer returns clinical chemistry alongside it. \code{"metab"}
#'   does not imply \code{"metab-t-clinical"}.
#' @param output_file a file path if desired, to autoatically save output in the standard size of height 1.54, width 1.225, recommended if output is known to be a single plot, not faceted
#' @param scale_factor simple way to scale plot if larger sizes are needed for posters or talks. recommend integers only
#' @param color_time_labels toggle TRUE/FALSE allows color tiles corresponding to time point colors to be used instead of x-axis
#' @param include_legend toggle TRUE/FALSE to include a legend
#' @param legend_position if include_legend == TRUE, position can be selected
#' @param verbose logical; toggle to include verbose information
#' @param epigen logical; toggle to include epigenetic features (only atac offered for this function). If you include a gene name, this could result in many many epigenetic features mapping to the one object.
#'
#' @returns a ggplot
#' @export plot_single_feature
#'
#' @author Daniel Katz and Christopher Jin
#'
#' @importFrom scales label_number
#' @importFrom dplyr mutate filter if_else c_across
#' @importFrom stringr str_to_sentence str_c
#' @importFrom ggplot2 ggplot
#'
#'
#' @examples
#' \dontrun{
#' plot_single_feature_test(feature="CAR 10:0",
#'   p_level = 0.05,
#'   selected_tissues = "all",
#'   selected_omes = "all",
#'   output_file = paste0("~/single_plot_test_scale.pdf"),
#'   scale_factor = 1,
#'   color_time_labels = F) +
#'theme(legend.position = "bottom")
#'
#'
#'plot_single_feature_test(feature="TFEB",
#'  p_level = 0.05,
#'  selected_tissues = "muscle",
#'  selected_omes = "prot-pr",
#'  output_file = paste0("~/single_plot_test_scale_tfeb_prot.png"),
#'  scale_factor = 1,
#'  color_time_labels = T,
#'  include_legend = F)
#'  }
plot_single_feature = function(feature,
                               selected_tissues = "all",
                               selected_omes = "all",
                               p_level = 0.05,
                               output_file = NULL,
                               scale_factor = 1,
                               color_time_labels = FALSE,
                               include_legend = TRUE,
                               legend_position = "right",
                               verbose = TRUE,
                               epigen = FALSE){

  selected_tissues = match.arg(selected_tissues,
                               choices = c("all", "blood", "adipose", "muscle"),
                               several.ok = TRUE)
  if(all(selected_tissues == 'all'))
    selected_tissues = MotrpacHumanPreSuspensionAnalysis::tissue_available_list()

  selected_omes <- match.arg(arg = tolower(selected_omes),
                             choices = c("all",
                                         "transcript-rna-seq",
                                         "prot-pr", "prot-ol",
                                         "prot-ph", "metab",
                                         "epigen-methylcap-seq",
                                         "epigen-atac-seq",
                                         MotrpacHumanPreSuspensionAnalysis::clinical_ome_list()),
                             several.ok = TRUE)

  if(all(selected_omes == 'all')) selected_omes = c("transcript-rna-seq",
                                                    "prot-pr", "prot-ol",
                                                    "prot-ph", "metab",
                                                    "epigen-methylcap-seq",
                                                    "epigen-atac-seq",
                                                    MotrpacHumanPreSuspensionAnalysis::clinical_ome_list())

  # Clinical chemistry is published as its own omes, one metabolomics and one
  # proteomics, and the loaders gate it behind `load_clinical`. It is loaded here
  # unconditionally and then filtered like every other ome, so `selected_omes` decides
  # whether it is plotted. `metab` does not imply `metab-t-clinical`, matching how
  # load_differential_analysis() exempts the clinical platform when it folds the metab
  # platforms together. This is the v2.0 split of clinical chemistry: the one
  # CLIN_CHEMISTRY_DA object became one per assay, and both arrive through the loaders.
  #
  # `metab` names a family rather than an assay: the loaders stack every research platform
  # under it, and the fold below replaces it with the platform in `assay`. The request is
  # expanded the same way so that the filtering happens in one vocabulary.
  # metab-t-clinical is not part of that family and is only ever matched by name.
  if ("metab" %in% selected_omes)
    selected_omes = unique(c(selected_omes,
                             MotrpacHumanPreSuspensionAnalysis::metab_only_list()))

  # the omes that actually measure this analyte, empty for anything that is not
  # clinical chemistry. Kept separate from the request so the error below can say
  # which ome would have to be asked for.
  feature_clinical_omes = .clinical_omes_for_feature(feature)
  plotted_clinical_omes = base::intersect(feature_clinical_omes, selected_omes)

  clin_chem_id = .is_clinical_chemistry_feature(feature)

  # Match feature to gene symbol and assay
  feature_info = MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
    #ok new as of 7/15/2025 -> now case insensitive.
    dplyr::filter(tolower(feature_id) == tolower(feature) |
                    tolower(gene_symbol) == tolower(feature) |
                    tolower(refmet_name) == tolower(feature)) %>%
    dplyr::semi_join(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE, by = "gene_symbol") %>%
    dplyr::mutate(feature_id = case_when(
      !is.na(refmet_name) ~ refmet_name,
      TRUE ~ feature_id
    ))

  gene_symbol_options <- as.character(unique(feature_info$gene_symbol))
  refmet_options <- as.character(unique(feature_info$refmet_name))
  label_options <- union(gene_symbol_options,refmet_options)
  label_options <- label_options[!is.na(label_options)]
  feature_label = label_options[1]
  if (is.na(feature_label) && !is.na(clin_chem_id)) feature_label = clin_chem_id

  if(length(label_options) > 1 & verbose) {
    message("if there are multiple gene symbols or refmet names corresponding
            to the input, the first is chosen for the purpose of labeling, use caution")
  }

  # Every non-epigenetic tissue and ome is loaded in one call and filtered afterwards.
  # The objects are lazily loaded and the load is cheap, and one unconditional load in
  # one vocabulary is far simpler than assembling the request ome by ome: clinical
  # chemistry stops being a special case appended after the fact, and honouring
  # `selected_omes` becomes an ordinary filter. Epigenomics stays behind `epigen`,
  # because it is fetched from AWS rather than lazily loaded.
  da_object = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
      selected_omes = "all",
      selected_tissues = "all",
      single_matrix = TRUE,
      epigen = epigen,
      load_clinical = TRUE,
      verbose = verbose
    ) %>%
    .fold_metab_platform_into_assay()

  if(verbose){
    message("DA is loaded automatically using requested settings.
              If any metab platform was requested, the metab features will default to
              any metabolomics platform measured.")
  }

  # Clinical chemistry analytes are not in the feature-to-gene table, so the clinical
  # id joins the ids matched there rather than replacing them.
  #
  # as.character() is load bearing: feature_id is a factor over every feature in the
  # study, and c() on a factor drops it to the integer codes behind its levels, so the
  # ids silently become row numbers that match nothing.
  requested_feature_ids = unique(c(as.character(feature_info$feature_id),
                                   clin_chem_id))
  requested_feature_ids = requested_feature_ids[!is.na(requested_feature_ids)]

  # Subset to relevant feature-specific data
  feature_specific_da = da_object %>%
    dplyr::filter(assay %in% selected_omes) %>%
    dplyr::filter(tissue %in% selected_tissues) %>%
    dplyr::filter(contrast_type == "exercise_with_controls") %>%
    dplyr::filter(feature_id %in% requested_feature_ids) %>%
    dplyr::filter(!is.na(Timepoint))

  if(nrow(feature_specific_da) == 0) {
    # a clinical analyte held back by the ome gate is a different problem from a
    # feature that is not in the data at all, and the fix is different too
    if(length(feature_clinical_omes) > 0 & length(plotted_clinical_omes) == 0) {
      stop(feature, " is measured by clinical chemistry (",
           paste(feature_clinical_omes, collapse = ", "),
           "), which is not among the omes you requested.
         Add it to `selected_omes`, or use \"all\", to plot it.")
    }
    stop("No differential analysis corresponds to your requested feature.
         Please double check your input matches something in the feature to gene
         mapping's feature_id column, or, for metabolites, the refmet_name column")
  }

  #assay name conversion table
  assay_names_table<-MotrpacHumanPreSuspensionAnalysis::assay_codes %>%
    dplyr::select(any_of(c('assay_code',
                           'assay_name',
                           'assay_abbreviation',
                           'omics_code',
                           'submission_code',
                           'ome_text',
                           'assay_short_text')
    )
    ) %>%
    distinct()

  # Loaded and filtered the same way as the differential analysis above, and filtered on
  # the same request: the join below is a full join, so an ome or a tissue left in here
  # would add rows for something that was not asked for.
  summary_stats = MotrpacHumanPreSuspensionAnalysis::load_summary_stats(selected_tissues = "all",
                                                                        selected_omes = "all",
                                                                        single_matrix = TRUE,
                                                                        load_clinical = TRUE,
                                                                        verbose = verbose) %>%
    .fold_metab_platform_into_assay() %>%
    dplyr::filter(assay %in% selected_omes) %>%
    dplyr::filter(tissue %in% selected_tissues) %>%
    dplyr::filter(feature_id %in% feature_specific_da$feature_id) %>%
    dplyr::mutate(SE = SD/sqrt(Count),
                  CI_95 = qt((1 + 0.95)/2, Count - 1))
  #this code is now matching the previous `mean_cl_normal` implementation, see: `Hmisc::smean.cl.normal`
  #where instead of using a strict wald CI, the SE multiplier is estimated from a t-distribution
  #makes the bounds slightly larger in most cases. Bigger diff with smaller n

  # Merge QC, pheno, and DA data for plotting. No longer have sample level results but
  local_feature_data = feature_specific_da %>%
    dplyr::full_join(summary_stats, by = c("tissue", "assay", "Timepoint", "randomGroupCode","feature_id")) %>%
    dplyr::mutate(
      below_p_cutoff = ifelse(adj_p_value >= p_level | is.na(adj_p_value), "above", "Below p threshold"),
      below_p_cutoff = factor(below_p_cutoff, levels = c("above", "Below p threshold"))
    ) %>%
    dplyr::filter(tissue %in% selected_tissues) %>%
    dplyr::left_join(assay_names_table,by = c("assay" = "assay_code")) %>%
    dplyr::mutate(
      # neither the clinical omes nor the conventional metabolomics platform has an
      # assay_codes row to join against, so the facet strip would read NA without
      # these. metab-t-conv is the same measurement as metab-t-clinical on a log2
      # scale, and the label says so to keep the two panels apart.
      assay_short_text = dplyr::if_else(
        assay %in% c("clinical-chemistry", "metab-t-clinical", "prot-clinical"),
        "Clin. Chem.", assay_short_text),
      assay_short_text = dplyr::if_else(
        assay == "metab-t-conv", "Conv. Metab (log2)", assay_short_text),
      tissue = stringr::str_to_sentence(tissue),
      Timepoint = dplyr::recode(Timepoint,
                                "pre_exercise" = "Pre",
                                "during_20_min" = "D20M",
                                "during_40_min" = "D40M",
                                "post_10_min" = "P10M",
                                "post_15_30_45_min" = "P15-45M",
                                "post_3.5_4_hr" = "P3.5/4H",
                                "post_24_hr" = "P24H"
      ),
      tissue_assay = stringr::str_c(tissue, " ", assay_short_text)
    )

  label_map <- c(
    "Below p threshold" = paste0("adj p < ", p_level),
    "above" = paste0("adj p >= ", p_level),
    "ADUControl" = "CON",
    "ADUEndur" = "EE",
    "ADUResist" = "RE"
  )

  # keyed on what is actually plotted, not on whether the feature has a clinical
  # measurement somewhere: a clinical analyte plotted from a research platform alone is
  # on the log2 scale like anything else. Only the clinical omes are absolute —
  # metab-t-conv is the conventional panel already log2 transformed.
  y_label = if (any(feature_specific_da$assay %in%
                    MotrpacHumanPreSuspensionAnalysis::clinical_ome_list())) {
    "Clinical chemistry features are presented in absolute scale. Others are log2(normalized value)"
  } else {
    "log2(normalized value)"
  }

  sc = scale_factor * 0.7
  # bounds = ""
  # `show.legend = TRUE` on every layer that contributes a key, together with the
  # fixed key sets on the two discrete scales below, makes the legend depend only on
  # the scale and not on which levels a given tissue happens to contain. Without it a
  # tissue with no significant timepoints draws a legend that differs from one that
  # has them, and `patchwork::plot_layout(guides = "collect")` only merges guides that
  # are identical, so it stacks the two variants instead of collecting them into one.
  g2 = ggplot(local_feature_data, aes(x = Timepoint, color = randomGroupCode))+
    geom_line(aes(y = Mean, group = randomGroupCode), linewidth = 0.5 * sc,
              show.legend = TRUE) +
    geom_point(aes(y = Mean, group = randomGroupCode),size = 1.8 * sc,
               show.legend = TRUE) +
    geom_errorbar(
      aes(
        #this code is now matching the previous `mean_cl_normal` implementation,
        #where instead of using a strict wald CI, the SE multiplier is estimated from a t-distribution
        #makes the bounds slightly larger in most cases. Bigger diff with smaller n
        ymin = Mean - CI_95*SE,
        ymax = Mean + CI_95*SE,
        group = randomGroupCode
      ),
      width = 0.2 * sc,
      linewidth = 0.4 * sc,
      alpha = 0.6,
      show.legend = TRUE
    ) +
    geom_point(
      aes(
        y = Mean,
        group = randomGroupCode,
        fill = below_p_cutoff
      ),
      size = 1.7 * sc,
      shape = 21,
      stroke = 0.1 * sc,
      color = "black",
      show.legend = TRUE
    ) +
    scale_color_manual(
      values = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS,
      labels = label_map,
      limits = sort(names(MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS))
    ) +
    scale_fill_manual(
      values = c(`Below p threshold` = "black", above = "white"),
      labels = label_map,
      drop = FALSE
    ) +
    facet_wrap(~ tissue_assay + feature_id, scales = "free_y") +
    ggtitle(feature_label) +
    xlab("Time point") +
    ylab(y_label) +
    theme_bw() +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5 * sc, color = "black"),
      axis.text.y = element_text(size = 7 * sc, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8 * sc),
      plot.title = element_text(size = 10 * sc, margin = margin(0, 0, 0, 0)),
      strip.text.x = element_text(size = 8 * sc, margin = margin(0.05 * sc, 0, 0.05 * sc, 0, "cm")),
      legend.text = element_text(size = 7 * sc, margin = margin(0, 0, 0, 0)),
      legend.spacing.x = unit(0.01 * sc, "in"),
      legend.spacing.y = unit(0.1 * sc, "in"),
      legend.box.spacing = unit(0.01 * sc, "in"),
      legend.key.size = unit(0.1 * sc, "in"),
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      axis.ticks = element_line(linewidth = 0.3 * sc),
      panel.grid = element_line(linewidth = 0.3 * sc),
      panel.border = element_rect(linewidth = 0.3 * sc),
      plot.margin = margin(0.05 * sc, 0.05 * sc, 0.05 * sc, 0.05 * sc, "in")
    )


  pb <- ggplot_build(g2)
  # Optional colored x-axis strip instead of x-axis labels
  if(color_time_labels == TRUE & length(pb$layout$layout$PANEL)==1) {
    layer_scales(g2)$y$range$range[1] -> g2_min
    layer_scales(g2)$y$range$range[2] -> g2_max
    abs(g2_max - g2_min) -> g2_range

    timepoint_colors <- c(
      "Pre" = "#BEBEBE",
      "D20M" = "#FDE725",
      "D40M" = "#BAD071",
      "P10M" = "#D1BBD7",
      "P15-45M" = "#AE76A3",
      "P3.5/4H" = "#882E72",
      "P24H" = "#61194F"
    )


    label_strip_df <- data.frame(
      Timepoint = factor(names(timepoint_colors), levels = names(timepoint_colors), ordered = TRUE),
      y = g2_min - 0.05 * g2_range
    ) %>%
      dplyr::filter(Timepoint %in% local_feature_data$Timepoint)

    g2 <- g2 +
      geom_tile(data = label_strip_df, aes(x = Timepoint, y = y, fill = Timepoint),
                height = 0.05 * g2_range, inherit.aes = FALSE) +
      scale_fill_manual(values = c(timepoint_colors, `Below p threshold` = "black", above = "white"),
                        labels = label_map) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  } else if (color_time_labels == TRUE & length(pb$layout$layout$PANEL)>1){
    if(verbose){
      message("Color strip will not render properly on plots with more than 1 facet.
            Suggest plotting each feature individually if color strip x-axis required ")
    }
  }

  #check if a legend is requested, and only include it if color-tile x-axis also not requested
  if(include_legend == TRUE & !is.null(legend_position) & color_time_labels == FALSE) {
    g2<-g2 + theme (legend.position = legend_position)
  }

  # Save the plot if output_file is provided
  if(length(unique(local_feature_data$selected_omes)) * length(unique(local_feature_data$tissue)) != 1 &
     is.null(output_file)) {
    if(verbose){
      message("Saved plot not requested or multi-faceted plot detected.
              Will not automatically save a plot. Please save plot manually.
              Recommended height per plot is 1.54 inches, and recommended width per plot
              of 1.225 inches with scale_factor set to 1.
              Use width of 1.715 if a right hand legend is desired")
    }
  } else {
    if(include_legend == TRUE & color_time_labels == FALSE) {
      #the plot needs to be wider if a legend is included
      ggsave(
        g2, filename = output_file,
        height = 2.2 * sc, width = 2.45 * sc, dpi = 600, units = "in"
      )
    } else {
      ggsave(
        g2, filename = output_file,
        height = 2.2 * sc, width = 1.75 * sc, dpi = 600, units = "in"
      )
    }
  }

  return(g2)

}


#' Name the metabolomics platform in the assay column of a summary statistics table
#'
#' The \code{*_DA} objects stack every metabolomics platform under \code{assay = "metab"}
#' and name the platform in its own column, and the pipeline now builds the
#' \code{*_SUM_STATS} objects the same way. Everything here keys on the platform, so it is
#' folded into \code{assay} and the redundant column dropped, which also keeps the join
#' against the differential analysis from meeting \code{platform} on both sides.
#'
#' Objects built before that change name the platform in \code{assay} already and carry no
#' \code{platform} column; they are returned untouched, so both shapes come out the same.
#'
#' @param x a data frame of summary statistics
#' @returns \code{x} with the platform named in \code{assay} and no \code{platform} column
#' @keywords internal
#' @noRd

.fold_metab_platform_into_assay = function(x) {
  if (!"platform" %in% colnames(x)) return(x)
  out = x %>%
    dplyr::mutate(assay = ifelse(assay == "metab",
                                 as.character(platform),
                                 as.character(assay))) %>%
    dplyr::select(-platform)
  return(out)
}


#' Look a feature up in the clinical chemistry summary statistics
#'
#' Case-insensitive lookup of the requested feature against the
#' \code{feature_id} columns of \code{BLOOD_METAB_T_CLINICAL_SUM_STATS} and \code{BLOOD_PROT_CLINICAL_SUM_STATS}.
#'
#' @param feature character; the feature name to look up
#' @returns A data frame of the matching \code{feature_id} and \code{assay} pairs,
#'   with zero rows when the feature is not a clinical chemistry analyte
#' @keywords internal
#' @noRd

.clinical_chemistry_matches = function(feature) {
  ids = dplyr::bind_rows(
    MotrpacHumanPreSuspensionAnalysis::BLOOD_METAB_T_CLINICAL_SUM_STATS,
    MotrpacHumanPreSuspensionAnalysis::BLOOD_PROT_CLINICAL_SUM_STATS
  ) %>%
    .fold_metab_platform_into_assay() %>%
    dplyr::select(dplyr::any_of(c("feature_id", "assay"))) %>%
    dplyr::distinct()

  return(ids[tolower(ids$feature_id) == tolower(feature), , drop = FALSE])
}


#' Check if a feature is a clinical chemistry analyte
#'
#' @param feature character; the feature name to look up
#' @returns The matched \code{feature_id} string (case-correct) if found,
#'   otherwise \code{NA_character_}
#' @keywords internal
#' @noRd

.is_clinical_chemistry_feature = function(feature) {
  matched = .clinical_chemistry_matches(feature)
  if (nrow(matched) > 0) return(matched$feature_id[1])
  return(NA_character_)
}


#' Which clinical chemistry omes carry a feature
#'
#' @param feature character; the feature name to look up
#' @returns A character vector of clinical omes measuring the feature, empty when
#'   the feature is not a clinical chemistry analyte
#' @keywords internal
#' @noRd

.clinical_omes_for_feature = function(feature) {
  return(unique(as.character(.clinical_chemistry_matches(feature)$assay)))
}
