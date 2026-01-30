#' @title Principal component analysis for the pre-cawg group
#' @description Basically just a wrapper function to plot with shapes and colors
#'    that fits the pre-cawg figure guidelines#'
#'
#' @param data_to_pca The actual numerical data that should be included in analysis.
#'    must be able to be converted into a numerical matrix
#' @param sample_metadata The correponding metadata to \code{data_to_pca}. The
#'    rownames of this matrix must correspond to colnames of the data matrix.
#' @param to_scale logical; Whether to scale the data when performing PCA.
#' @param custom_title character; whatever title is desired for the figure
#'
#' @return a ggplot object with the desired PCA plotting aesthetic.
#' @export plot_precovid_pca
#'

plot_precovid_pca = function(data_to_pca,
                             sample_metadata,
                             to_scale = TRUE,
                             custom_title = ""){
  if(length(intersect(rownames(sample_metadata), colnames(data_to_pca))) == 0)
    stop("You may have forgotten to specify rownames of the metadata that match
         the columns of the data")
  sample_metadata = sample_metadata[match(colnames(data_to_pca), rownames(sample_metadata)), ] #reorder so they're in the same order as colnames of raw_counts
  pca = prcomp(t(data_to_pca), scale. = to_scale, center = T)
  pca_coor = pca$x[,1:2] %>% as.data.frame()
  pca_coor$vialLabel = rownames(pca_coor)
  var_exp = summary(pca)[["importance"]][2,c(1,2)]*100
  pca_metadata = merge(sample_metadata, pca_coor, by = "vialLabel")
  #note: we don't have shapes and stuff already determined so I'm just letting it go default mode for shape
  pca_plot = pca_metadata %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(
      aes(color = Timepoint, shape = randomGroupCode),
      size = 2,
      stroke = 0.5,
      alpha = 0.9
    ) +
    scale_color_manual(values = MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS) +
    scale_shape_manual(
      values = c("ADUEndur" = 21,
                 "ADUControl" = 23,
                 "ADUResist" = 25)
    ) +
    labs(
      title = custom_title,
      x = paste0('PC1 (', var_exp[1], '%)'),
      y = paste0('PC2 (', var_exp[2], '%)')
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )

  return(pca_plot)

}
