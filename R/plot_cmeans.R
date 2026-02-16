#' @title Save a Plot of Fuzzy C-Means Cluster Trajectories
#'
#' @description Essentially, a ggplot2 version of \code{Mfuzz::plot.mfuzz2} with
#'   a different color scheme that saves a plot to a file, rather than
#'   displaying it.
#'
#' @param FCM one of the elements of the list produced by
#'   \code{\link{run_cmeans}}. See examples.
#' @param keep_clusters integer; the identifier numbers of clusters to include
#'   in the plot. Defaults to all clusters.
#' @param min_membership numeric; to be included in the plot, the probability
#'   that a feature belongs to a cluster should be \eqn{\geq}
#'   \code{min_membership}.
#' @param common_ylim logical; whether the plots for each cluster should have
#'   the same y-axis limits.
#' @param ncol integer; the maximum number of clusters to place on a single row.
#'   The number of rows is determined automatically.
#' @inheritParams ggplot2::ggsave
#' @param modality character; which exercise modalities should be used for
#'   plotting? One of \code{"both"}, \code{"Endur"}, or \code{"Resist"}.
#'
#' @author Tyler Sagendorf
#'
#' @seealso \code{\link{run_cmeans}}
#'
#' @import ggplot2
#' @importFrom dplyr %>% mutate across inner_join arrange filter
#' @importFrom ggpubr ggarrange annotate_figure
#' @importFrom grid textGrob gpar
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer nest
#'
#' @export plot_cmeans
#'
#' @examples
#' \dontrun{
#' # FCM results
#' FCM <- run_cmeans(selected_tissues = "adipose")
#'
#' plot_cmeans(FCM = FCM[["adipose"]],
#'             filename = "adipose_clusters.pdf",
#'             min_membership = 0.3)
#' }

plot_cmeans <- function(FCM,
                        filename,
                        keep_clusters,
                        min_membership = 0,
                        common_ylim = TRUE,
                        ncol = 4L,
                        dpi = 250,
                        modality = c("both", "Endur", "Resist"))
{
  if (!inherits(FCM, "fclust")) {
    stop("`FCM` must be an object of class 'fclust'. See documentation ",
         "for details.")
  }

  min_membership <- max(0, min(min_membership, 1, na.rm = TRUE))

  modality <- match.arg(modality,
                        choices = c("both", "Endur", "Resist"))

  if (modality == "both")
    modality <- c("Endur", "Resist")

  if (missing(keep_clusters)) {
    keep_clusters <- as.integer(rownames(FCM[["centers"]]))
  }

  mem <- FCM[["membership"]][, keep_clusters, drop = FALSE]
  centers <- FCM[["centers"]][keep_clusters, , drop = FALSE]
  colnames(mem) <- rownames(centers) <- paste0("Cluster ", colnames(mem))

  colnames(centers) <- sub("group_timepointADU([^ ]+) .*", "\\1",
                           colnames(centers))
  colnames(centers) <- paste0("timepoint_", colnames(centers))

  zmat <- FCM[["input"]]
  colnames(zmat) <- sub("group_timepointADU([^ ]+) .*", "\\1", colnames(zmat))
  colnames(zmat) <- paste0("timepoint_", colnames(zmat))

  if (length(modality) == 1L) {
    centers <- centers[, grepl(modality, colnames(centers))]
    zmat <- zmat[, grepl(modality, colnames(zmat))]
  }

  modalities <- c(
    c("Endur")[any(grepl("Endur", colnames(centers)))],
    c("Resist")[any(grepl("Resist", colnames(centers)))]
  )

  ## Include pre_exercise timepoints for both exercise modalities
  pre_ex <- matrix(data = 0, nrow = nrow(zmat), ncol = length(modality))
  rownames(pre_ex) <- rownames(zmat)
  colnames(pre_ex) <- sprintf("timepoint_%s.pre_exercise",
                              modality)

  zmat <- cbind(pre_ex, zmat)
  zmat <- zmat[, c(grep("Endur", colnames(zmat)),
                   grep("Resist", colnames(zmat)))]

  # Include two columns of zeros in the centers
  pre_ex_center <- matrix(data = 0, nrow = nrow(centers),
                          ncol = length(modalities))
  rownames(pre_ex_center) <- rownames(centers)
  colnames(pre_ex_center) <- sprintf("timepoint_%s.pre_exercise",
                                     modalities)

  centers <- cbind(centers, pre_ex_center)
  centers <- centers[, c(grep("Endur", colnames(centers)),
                         grep("Resist", colnames(centers)))]

  cluster_assignment <- data.frame(feature = names(FCM$cluster),
                                   cluster = paste0("Cluster ", FCM$cluster))

  z_df <- zmat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(cols = -feature,
                        names_to = "timepoint",
                        names_pattern = "timepoint_(.*)",
                        values_to = "z") %>%
    dplyr::mutate(modality = sub("\\..*", "", timepoint),
           timepoint = sub("^[^.]+\\.(.*)$", "\\1", timepoint),
           timepoint = gsub("(?<=\\d)_(?=\\d)", "/", timepoint, perl = TRUE),
           timepoint = gsub("_", " ", timepoint),
           across(.cols = c(modality, timepoint),
                  .fns = ~ factor(.x, levels = unique(.x))))

  cluster_df <- mem %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(cols = -feature,
                        names_to = "cluster",
                        values_to = "membership") %>%
    dplyr::filter(membership >= min_membership) %>%
    dplyr::inner_join(cluster_assignment, by = c("feature", "cluster")) %>%
    dplyr::mutate(cluster = factor(cluster, levels = rownames(centers)))

  df <- inner_join(z_df, cluster_df, by = "feature")  %>%
    dplyr::arrange(cluster, membership) %>%
    dplyr::mutate(feature = factor(feature, levels = unique(feature))) %>%
    droplevels.data.frame() %>%
    dplyr::mutate(cluster = factor(cluster, levels = rownames(centers)))

  # Cluster centroids
  center_df <- centers %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cluster") %>%
    tidyr::pivot_longer(cols = -cluster,
                        names_to = "timepoint",
                        names_pattern = "timepoint_(.*)",
                        values_to = "center") %>%
    dplyr::mutate(modality = sub("\\..*", "", timepoint),
           modality = factor(modality, levels = levels(z_df$modality)),
           timepoint = sub("^[^.]+\\.(.*)$", "\\1", timepoint),
           timepoint = gsub("(?<=\\d)_(?=\\d)", "/", timepoint, perl = TRUE),
           timepoint = gsub("_", " ", timepoint),
           timepoint = factor(timepoint, levels = unique(timepoint)),
           cluster = factor(cluster, levels = rownames(centers)))

  if (common_ylim) {
    ylim <- range(df$z)
  } else {
    ylim <- rep(NA, 2L)
  }

  df_list <- tidyr::nest(df, .by = cluster)
  center_df_list <- tidyr::nest(center_df, .by = cluster)

  ## Iridescent color scheme from https://personal.sron.nl/~pault/
  cmeans_palette <- c("#fefbe9", "#fcf7d5", "#f5f3c1", "#eaf0b5", "#ddecbf",
                      "#d0e7ca", "#c2e3d2", "#b5ddd8", "#a8d8dc", "#9bd2e1",
                      "#8dcbe4", "#81c4e7", "#7bbce7", "#7eb2e4", "#88a5dd",
                      "#9398d2", "#9b8ac4", "#9d7db2", "#9a709e", "#906388",
                      "#805770", "#684957", "#46353a")

  ## List of plots for each cluster
  plotlist <- lapply(seq_len(nrow(df_list)), function(i) {
    df_i <- df_list$data[[i]]
    centroid_i <- center_df_list$data[[i]]

    pi <- ggplot(df_i, aes(x = timepoint, y = z)) +
      geom_line(aes(color = membership, group = feature)) +
      geom_line(aes(x = timepoint, y = center, group = 1),
                data = centroid_i, color = "black") +
      scale_x_discrete(name = NULL,
                       expand = expansion(add = 0.2)) +
      scale_y_continuous(limits = ylim) +
      scale_color_gradientn(name = "Membership\nProbability",
                            colors = cmeans_palette,
                            values = seq(0, 1,
                                         length.out = length(cmeans_palette)),
                            limits = c(0, 1),
                            breaks = seq(0, 1, 0.2)) +
      labs(y = NULL,
           title = df_list$cluster[i]) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            strip.text = element_text(face = "bold",
                                      margin = margin(t = 3)),
            plot.title = element_text(size = rel(0.9),
                                      face = "bold", hjust = 0.5,
                                      margin = margin(b = 0)),
            strip.background = element_blank(),
            axis.title.y = element_text(color = "black", face = "bold",
                                        margin = margin(r = 8)),
            axis.text.x = element_text(size = rel(0.7),
                                       color = "black", angle = 90,
                                       vjust = 0.5, hjust = 1),
            axis.text.y = element_text(color = "black"))

    if (length(modality) == 2L)
      pi <- pi +
        facet_grid(~ modality)

    return(pi)
  })

  ncol <- min(length(plotlist), ncol)
  nrow <- ceiling(length(plotlist) / ncol)
  # Each cluster occupies a space that is 2x2 inches
  height <- nrow * 2
  width <- ncol * (1.5 + 0.5 * (length(modality) == 2L)) + 0.5 # 0.5 = space for legend

  p <- ggarrange(plotlist = plotlist,
                 ncol = ncol,
                 nrow = nrow,
                 common.legend = TRUE,
                 legend = "right")

  # Add global y-axis title
  p <- annotate_figure(
    p = p,
    left = textGrob("Scaled Z-Score", rot = 90, vjust = 1,
                    gp = gpar(fontsize = 11,
                              fontface = "bold"))
  )

  # Save to file
  ggsave(filename = filename, plot = p,
         height = height, width = width,
         dpi = dpi, bg = "white")
}
