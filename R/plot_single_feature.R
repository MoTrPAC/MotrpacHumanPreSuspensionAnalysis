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
#' @param repo_local_dir Deprecated and ignored; passed through to
#'   \code{\link{load_differential_analysis}}, which prints a message if it is
#'   supplied.
#' @param verbose logical; toggle to include verbose information
#' @param epigen logical; toggle to include epigenetic features (only atac offered for this function). If you include a gene name, this could result in many many epigenetic features mapping to the one object. Only the epigenomic files for the requested tissues and omes are downloaded.
#' @param qc_data optional; the nested list returned by
#'   \code{MotrpacHumanPreSuspensionData::load_qc()}, for users with access to the
#'   individual-level data. The shipped \code{*_SUM_STATS} objects cover only the
#'   features in the differential analysis, and for the epigenomic omes only the
#'   features with \code{adj_p_value < 0.05} in at least one contrast. When a
#'   requested feature has differential analysis but no summary statistics, they are
#'   computed from \code{qc_data} the way the shipped objects are: \code{qc_norm}
#'   restricted to \code{visitcode == "ADU_BAS"} samples, summarised per
#'   \code{randomGroupCode} and \code{Timepoint}. Load it with the tissues and omes
#'   being plotted, and \code{epigen = TRUE} for the epigenomic omes.
#'
#' @section Features without summary statistics:
#' The lines, points and error bars are drawn from the summary statistics, not from
#' the differential analysis. A feature that is in the differential analysis but not in
#' the summary statistics, most often an epigenomic feature that is not significant in
#' any contrast, has nothing to draw. The function then says so with a
#' \code{message()} naming the features and returns a plot with empty panels for them,
#' unless \code{qc_data} is supplied to compute the missing statistics.
#'
#' @returns a ggplot
#' @export plot_single_feature
#'
#' @author Daniel Katz and Christopher Jin
#'
#' @importFrom scales label_number
#' @importFrom dplyr mutate filter if_else c_across
#' @importFrom stringr str_to_sentence str_c
#' @import ggplot2
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
                               repo_local_dir = NULL,
                               verbose = TRUE,
                               epigen = FALSE,
                               qc_data = NULL){

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
  # `selected_omes` becomes an ordinary filter. Epigenomics is downloaded from the CDN,
  # one file per tissue and ome, so it is loaded separately and only for what was
  # requested.
  # named rather than "all": "all" with epigen = FALSE reports epigenomics as skipped
  da_object = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
      selected_omes = c("transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol", "metab",
                        MotrpacHumanPreSuspensionAnalysis::clinical_ome_list()),
      selected_tissues = "all",
      single_matrix = TRUE,
      epigen = FALSE,
      repo_local_dir = repo_local_dir,
      load_clinical = TRUE,
      verbose = verbose
    )

  epigen_omes = base::intersect(selected_omes, c("epigen-atac-seq", "epigen-methylcap-seq"))
  # only the tissue and ome pairs that were measured have a file to download
  epigen_measured = MotrpacHumanPreSuspensionAnalysis::OME_TISSUE_CODE %>%
    dplyr::filter(ome %in% epigen_omes, tissue %in% selected_tissues)
  if (epigen && nrow(epigen_measured) > 0) {
    epigen_da = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
      selected_omes = unique(as.character(epigen_measured$ome)),
      selected_tissues = unique(as.character(epigen_measured$tissue)),
      single_matrix = TRUE,
      epigen = TRUE,
      verbose = verbose
    )
    da_object = data.table::rbindlist(list(da_object, epigen_da), use.names = TRUE, fill = TRUE)
    data.table::setorderv(da_object, cols = "contrast", order = 1L)
  } else if (!epigen && length(epigen_omes) > 0 && verbose) {
    message("You've requested one or more epigenetic omes (via explicit selection ",
            "or \"all\") but `epigen = FALSE`, so epigenetic data will be skipped. ",
            "Set `epigen = TRUE` to load epigenetic data.")
  }

  da_object = .fold_metab_platform_into_assay(da_object)

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
  assay_names_table<-MotrpacBicQC::assay_codes %>%
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
    as.data.frame()

  # Epigenomic summary statistics ship only for features significant in some contrast,
  # so a feature can have differential analysis and nothing to draw.
  plot_keys = c("tissue", "assay", "feature_id")
  missing_sum_stats = feature_specific_da %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(plot_keys))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
    dplyr::anti_join(summary_stats %>%
                       dplyr::distinct(dplyr::across(dplyr::all_of(plot_keys))) %>%
                       dplyr::mutate(dplyr::across(dplyr::everything(), as.character)),
                     by = plot_keys)

  if (nrow(missing_sum_stats) > 0 && !is.null(qc_data)) {
    # Timepoint is a factor on both tiers and its levels set the x-axis order; a
    # character column here would turn the joined column character and sort it
    # alphabetically.
    timepoint_levels = levels(summary_stats$Timepoint)
    summary_stats = dplyr::bind_rows(
      summary_stats %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(c(plot_keys, "Timepoint", "randomGroupCode")),
                                    as.character)),
      .summary_stats_from_qc(qc_data, missing_sum_stats,
                             timepoints = timepoint_levels, verbose = verbose)
    ) %>%
      dplyr::mutate(Timepoint = factor(Timepoint, levels = timepoint_levels))
    missing_sum_stats = dplyr::anti_join(missing_sum_stats,
                                         dplyr::distinct(summary_stats[, plot_keys]),
                                         by = plot_keys)
  }

  if (nrow(missing_sum_stats) > 0) {
    message("No summary statistics for ",
            paste0(missing_sum_stats$tissue, " ", missing_sum_stats$assay, " ",
                   missing_sum_stats$feature_id, collapse = ", "),
            ", so their panels are empty.
         Epigenomic summary statistics ship only for features with adj_p_value < 0.05
         in at least one contrast.",
            if (is.null(qc_data))
              " With individual-level data access, pass `qc_data` from
         MotrpacHumanPreSuspensionData::load_qc() to compute them."
            else
              " The `qc_data` you supplied does not carry these features or their acute samples.")
  }

  summary_stats = summary_stats %>%
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
      # metab-t-conv carries the same measurement as metab-t-clinical on a log2
      # scale. Upstream labels both from the LAB assay family, so the override keeps
      # the two panels distinguishable.
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


#' Name the metabolomics platform in the assay column
#'
#' Both the \code{*_DA} and the \code{*_SUM_STATS} objects stack every metabolomics
#' platform under \code{assay = "metab"} and name the platform in its own column.
#' Everything here keys on the platform, so it is folded into \code{assay} and the
#' redundant column dropped, which also keeps the join between the two tiers from
#' meeting \code{platform} on both sides.
#'
#' The \code{platform} column is required. Summary statistics built before v2.0.1 named
#' the platform in \code{assay} and carried no such column; on those objects every metab
#' row would silently keep the family name and match no differential analysis row, so
#' they are rejected here rather than plotted as a panel with no points.
#'
#' @param x a data frame of differential analysis results or summary statistics,
#'   combined across omes so that the metabolomics \code{platform} column is present
#' @returns \code{x} with the platform named in \code{assay} and no \code{platform} column
#' @keywords internal
#' @noRd

.fold_metab_platform_into_assay = function(x) {
  if (!"platform" %in% colnames(x)) {
    stop("These data objects carry no `platform` column, so they predate v2.0.1 of the
         package. Reinstall MotrpacHumanPreSuspensionAnalysis to get summary statistics
         that name the metabolomics platform the way the differential analysis does.")
  }
  out = x %>%
    dplyr::mutate(assay = ifelse(assay == "metab",
                                 as.character(platform),
                                 as.character(assay))) %>%
    dplyr::select(-platform)
  return(out)
}


#' Summary statistics for chosen features from load_qc() output
#'
#' Computes what the shipped \code{*_SUM_STATS} objects hold, for features they do not
#' carry: \code{qc_norm} restricted to \code{visitcode == "ADU_BAS"} samples, then
#' \code{Count}, \code{Mean} and \code{SD} per \code{randomGroupCode},
#' \code{Timepoint} and feature.
#'
#' @param qc_data the nested list returned by
#'   \code{MotrpacHumanPreSuspensionData::load_qc()}, tissue then ome
#' @param features a data frame with \code{tissue}, \code{assay} and
#'   \code{feature_id} columns; metabolomics rows name the platform in \code{assay}
#' @param timepoints character; the \code{Timepoint} values the summary statistics
#'   use. Samples at any other value, or at none, are left out.
#' @param verbose logical; report samples left out for their group or timepoint
#' @returns A data frame with \code{tissue}, \code{assay}, \code{randomGroupCode},
#'   \code{Timepoint}, \code{feature_id}, \code{Count}, \code{Mean} and \code{SD},
#'   with zero rows when \code{qc_data} carries none of the features
#' @keywords internal
#' @noRd

.summary_stats_from_qc = function(qc_data, features, timepoints = NULL, verbose = TRUE) {
  empty = data.frame(tissue = character(0), assay = character(0),
                     randomGroupCode = character(0), Timepoint = character(0),
                     feature_id = character(0), Count = integer(0),
                     Mean = numeric(0), SD = numeric(0))

  # load_qc() nests tissue, then ome, then qc_norm / sample_metadata. One level of it
  # passed on its own, e.g. qc$muscle, would otherwise match nothing without saying so.
  is_qc_entry = function(x) {
    return(is.list(x) && all(c("qc_norm", "sample_metadata") %in% names(x)))
  }
  is_tissue_list = function(tissue_list) {
    return(is.list(tissue_list) && length(tissue_list) > 0 &&
             all(vapply(tissue_list, is_qc_entry, logical(1))))
  }
  if (!is.list(qc_data) || is.null(names(qc_data)) || is_qc_entry(qc_data) ||
      !all(tolower(names(qc_data)) %in% c("adipose", "blood", "muscle")) ||
      !all(vapply(qc_data, is_tissue_list, logical(1)))) {
    stop("`qc_data` must be the whole nested list returned by
         MotrpacHumanPreSuspensionData::load_qc(): tissue, then ome, then qc_norm and
         sample_metadata. Pass the load_qc() result itself, not one tissue or ome of it.")
  }
  names(qc_data) = tolower(names(qc_data))

  out = list()
  for (combo in split(features, list(features$tissue, features$assay), drop = TRUE)) {
    tissue = combo$tissue[1]
    assay = combo$assay[1]
    qc_entry = qc_data[[tolower(tissue)]][[assay]]
    if (is.null(qc_entry)) next

    required = c("vialLabel", "visitcode", "randomGroupCode", "Timepoint")
    absent = setdiff(required, colnames(qc_entry$sample_metadata))
    if (length(absent) > 0) {
      stop("`qc_data` ", tissue, " ", assay, " sample_metadata has no ",
           paste(absent, collapse = ", "), " column.")
    }

    qc_norm = qc_entry$qc_norm
    qc_norm = qc_norm[rownames(qc_norm) %in% combo$feature_id, , drop = FALSE]
    # vialLabel is read in as numeric, and a numeric index selects columns by position.
    # distinct(): a vial listed twice would otherwise be counted twice.
    sample_metadata = qc_entry$sample_metadata %>%
      dplyr::mutate(vialLabel = as.character(vialLabel)) %>%
      dplyr::filter(visitcode == "ADU_BAS", vialLabel %in% colnames(qc_norm)) %>%
      dplyr::distinct(vialLabel, .keep_all = TRUE)
    unplaceable = is.na(sample_metadata$randomGroupCode) | is.na(sample_metadata$Timepoint) |
      (!is.null(timepoints) & !as.character(sample_metadata$Timepoint) %in% timepoints)
    if (any(unplaceable) && verbose) {
      message(sum(unplaceable), " acute ", tissue, " ", assay, " sample(s) in `qc_data` ",
              "have no group or a timepoint outside the summary statistics and are left out.")
    }
    sample_metadata = sample_metadata[!unplaceable, , drop = FALSE]
    qc_norm = qc_norm[, sample_metadata$vialLabel, drop = FALSE]
    if (nrow(qc_norm) == 0 || ncol(qc_norm) == 0) next
    if (!all(vapply(qc_norm, is.numeric, logical(1)))) {
      stop("`qc_data` ", tissue, " ", assay, " qc_norm is not numeric.")
    }

    out[[paste(tissue, assay)]] = data.frame(
      feature_id = rep(rownames(qc_norm), times = ncol(qc_norm)),
      randomGroupCode = rep(as.character(sample_metadata$randomGroupCode),
                            each = nrow(qc_norm)),
      Timepoint = rep(as.character(sample_metadata$Timepoint), each = nrow(qc_norm)),
      Value = unlist(qc_norm, use.names = FALSE)
    ) %>%
      dplyr::filter(!is.na(Value)) %>%
      dplyr::group_by(randomGroupCode, Timepoint, feature_id) %>%
      dplyr::summarize(Count = dplyr::n(), Mean = mean(Value), SD = stats::sd(Value),
                       .groups = "drop") %>%
      dplyr::mutate(tissue = tissue, assay = assay, .before = 1)
  }

  if (length(out) == 0) return(empty)
  return(as.data.frame(dplyr::bind_rows(out)))
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
