#' @title canCor analysis cor correlation to principal components
#'
#' @description Performs PCA and correlates PCs to selected metadata using
#' \code{variancePartition::canCorPairs}.
#'
#' @param data_to_pca The actual numerical data that should be included in analysis.
#'    must be able to be converted into a numerical matrix
#' @param sample_metadata The correponding metadata to \code{data_to_pca}. The
#'    rownames of this matrix must correspond to colnames of the data matrix.
#' @param to_scale logical; Whether to scale the data when performing PCA.
#' @param custom_title character; whatever title is desired for the figure
#' @param interesting_variables character; a vector of variables that are
#'    going to be correlated using canCors to the PCs. The # of PCs in
#'    \code{num_pcs} are included by default.
#' @param num_pcs numeric; the number of principal components desired to be included.
#'
#' @return a pheatmap object with the desired correlations to principal components
#' @export plot_precovid_cca
#' @author Christopher Jin
#'

plot_precovid_cca = function(data_to_pca,
                             sample_metadata,
                             to_scale = TRUE,
                             custom_title = "",
                             interesting_variables = c("pid", "randomGroupCode", "Timepoint", "Sex", "BMI", "calculatedAge"),
                             num_pcs = 5
){
  check_package_installation(pkg = "RColorBrewer")
  check_package_installation(pkg = "viridis")
  check_package_installation(pkg = "pheatmap")
  check_package_installation(pkg = "variancePartition")
  check_package_installation(pkg = "MotrpacHumanPreSuspensionData")

  if(length(intersect(rownames(sample_metadata), colnames(data_to_pca))) == 0)
    stop("You may have forgotten to specify rownames of the metadata that match
         the columns of the data")

  sample_metadata = sample_metadata[match(colnames(data_to_pca),
                                          rownames(sample_metadata)), ]
  pca = prcomp(t(data_to_pca), scale. = to_scale, center = TRUE)
  pca_coor = pca$x[, 1:num_pcs] %>% as.data.frame()
  pca_coor$vialLabel = rownames(pca_coor)
  var_exp = summary(pca)[["importance"]][2, c(1:num_pcs)] * 100
  pca_metadata = merge(sample_metadata, pca_coor, by = "vialLabel")

  interesting_variables = c(
    paste0("PC", 1:num_pcs),
    interesting_variables #from user inputted info
  )
  form <- reformulate(interesting_variables, response = NULL)
  c_pca = suppressWarnings(canCorPairs(form, pca_metadata))

  rownames(c_pca) = colnames(c_pca) = interesting_variables
  coul <- brewer.pal(9, "Reds")
  pcs = paste0("PC", 1:num_pcs)
  pc_labels = paste0(pcs, " (", sprintf("%.1f", var_exp), "%)")
  colnames(c_pca)[match(pcs, colnames(c_pca))] = pc_labels

  breaksList = seq(0, 1.00, by = 0.1)

  plot = pheatmap::pheatmap(c_pca[!rownames(c_pca) %in% pcs, pc_labels],
                            color = coul,
                            breaks = breaksList,
                            fontsize = 14,
                            cluster_cols = FALSE,
                            main = custom_title)

  return(plot)
}
