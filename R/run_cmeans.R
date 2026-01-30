#' @title Fuzzy C-Means (FCM) Clustering
#'
#' @inheritParams run_cameraPR
#' @param selected_omes character; one or more of \code{"transcript-rna-seq"},
#'   \code{"prot-pr"}, \code{"prot-ol"}, \code{"prot-ph"}, or \code{"metab"}.
#' @param num_clusters_adipose,num_clusters_blood,num_clusters_muscle integer;
#'   the number of clusters desired for each tissue. Defaults to 13 for adipose,
#'   13 for blood, and 12 for muscle. If more than one number is provided,
#'   \code{Mfuzz::Dmin} will be used to generate a plot of the minimum centroid
#'   distances, and the user will be asked to specify the optimal cluster number
#'   in the console for each tissue.
#' @param modality character; which exercise modalities should be used for FCM?
#'   One of \code{"both"}, \code{"Endur"}, or \code{"Resist"}.
#'
#' @returns A named list of objects where names are tissues. Each object is of
#'   class \code{"fclust"} with additional list component \code{"input"} for the
#'   matrix of scaled z-scores used as input for FCM (both modalities included,
#'   even when \code{modality != "both"}).
#'
#' @author Tyler Sagendorf, Christopher Jin
#'
#' @seealso \code{\link{plot_cmeans}}, \code{\link{run_cluster_cameraPR}},
#'   \code{\link{run_cluster_ORA}}
#'
#' @importFrom data.table rbindlist
#' @importFrom stats complete.cases
#' @importFrom Biobase ExpressionSet
#' @importFrom Mfuzz mestimate mfuzz Dmin
#'
#' @export run_cmeans
#'
#' @examples
#' \dontrun{
#'   x1 <- run_cmeans()
#'   names(x1) # list available components
#'
#'   # Try a range of cluster numbers for one tissue
#'   x2 <- run_cmeans(selected_tissues = "adipose",
#'                    num_clusters_adipose = 3:14)
#'
#'   # FCM for a single modality
#'   x3 <- run_cmeans(selected_tissues = "adipose",
#'                    modality = "Endur")
#' }

run_cmeans <- function(selected_tissues = c("all", "adipose",
                                            "blood", "muscle"),
                       selected_omes = c("transcript-rna-seq",
                                         "prot-pr", "prot-ol",
                                         "prot-ph", "metab"),
                       num_clusters_adipose = 13L,
                       num_clusters_blood = 13L,
                       num_clusters_muscle = 12L,
                       modality = c("both", "Endur", "Resist"))
{
  on.exit(gc())
  set.seed(0)
  check_package_installation(pkg = "Mfuzz", fun = "run_cmeans")
  check_package_installation(pkg = "Biobase", fun = "run_cmeans")

  modality <- match.arg(modality,
                        choices = c("both", "Endur", "Resist"))

  selected_omes <- match.arg(arg = tolower(selected_omes),
                             choices = c("transcript-rna-seq",
                                         "prot-pr", "prot-ol",
                                         "prot-ph", "metab"),
                             several.ok = TRUE)

  num_clusters <- list("adipose" = num_clusters_adipose,
                       "blood" = num_clusters_blood,
                       "muscle" = num_clusters_muscle)

  num_clusters <- lapply(num_clusters, function(nc) {
    if (!is.vector(nc, mode = "numeric")) {
      stop("One or more of `num_clusters_adipose`, `num_clusters_blood`, or ",
           "`num_clusters_muscle` are not integer vectors.")
    }

    sort(unique(as.integer(pmax(2L, nc, na.rm = TRUE))))
  })

  ## Prepare DA results ----
  DA_list <- .prepare_DA_results(DA_list = NULL,
                                 selected_omes = selected_omes,
                                 selected_tissues = selected_tissues,
                                 convert_features = FALSE,
                                 .contrast_type = "exercise_with_controls")

  # Remove adipose global and phosphoproteomics due to insufficient timepoints
  DA_list[c("adipose.prot-pr", "adipose.prot-ph")] <- NULL
  nm <- names(DA_list)
  DA_list <- lapply(names(DA_list), function(name_i) {
    zi <- DA_list[[name_i]]

    zi <- zi[, !grepl("during", colnames(zi))]
    ome_i <- sub("^[^.]+\\.", "", name_i)
    rownames(zi) <- paste(ome_i, rownames(zi))

    return(zi)
  })
  names(DA_list) <- nm

  tissue <- sub("\\..*", "", names(DA_list))

  # Nest by tissue and stack matrices
  DA_list <- split(do.call(list, DA_list), tissue)

  zmat_list <- lapply(DA_list, function(li) {
    zi <- do.call(what = rbind, args = li)
    zi <- zi[complete.cases(zi), ]

    return(zi)
  })

  # Create ExpressionSet objects for mfuzz. Each row is scaled by dividing by
  # the standard deviation calculated after including two columns of zeros to
  # represent the pre-exercise timepoints. If this is not done, features that
  # change similarly in all contrasts will move further away from 0.
  eset_list <- lapply(zmat_list, function(zi) {
    zero_mat <- matrix(data = 0, nrow = nrow(zi), ncol = 2L)
    colnames(zero_mat) <- paste0(c("Endur", "Resist"),
                                 ".pre_exercise")
    tmp <- cbind(zero_mat, zi)

    # Scale features
    sd <- apply(tmp, 1L, sd)

    zi <- sweep(zi, 1L, sd, FUN = "/")

    zi <- zi[order(sd, decreasing = TRUE), ]

    Biobase::ExpressionSet(assayData = zi)
  })

  cmeans_list <- vector(mode = "list", length = length(eset_list))
  names(cmeans_list) <- names(eset_list)

  for (tissue_i in names(eset_list)) {
    eset_i <- eset_list[[tissue_i]]

    # m is determined before subsetting, since it becomes too large with few
    # contrast columns.
    m_i <- Mfuzz::mestimate(eset_i)

    # Optionally subset to specific modality
    if (modality != "both") {
      eset_i <- eset_i[, grepl(modality, colnames(eset_i)), drop = FALSE]
    } else {
      eset_i <- eset_i[, c(grep("Endur", colnames(eset_i)),
                           grep("Resist", colnames(eset_i)))]
    }

    num_clusters_i <- num_clusters[[tissue_i]]

    if (length(num_clusters_i) > 1L) {
      message("Determining optimal number of clusters for ", tissue_i, "...")
      min_centroid_dist <- Mfuzz::Dmin(
        eset = eset_i,
        m = m_i,
        crange = num_clusters_i,
        repeats = 1L,
        visu = FALSE
      )

      plot(x = num_clusters_i, y = min_centroid_dist,
           xlab = "Number of clusters",
           ylab = "Min. centroid dist.",
           main = tissue_i)

      num_clusters_i <- ""

      while (num_clusters_i %in% c("", "0", "1")) {
        num_clusters_i <- readline(
          prompt = sprintf(
            "Enter the optimal number of clusters (>= 2) for %s: ",
            tissue_i
          )
        )
        num_clusters_i <- gsub("[^[:digit:].]", "", num_clusters_i)
        num_clusters_i <- sub("\\..*", "", num_clusters_i)
        num_clusters_i <- sub("^[0]+", "", num_clusters_i)
      }

      num_clusters_i <- as.integer(num_clusters_i)
    }

    # FCM clustering
    FCM_i <- Mfuzz::mfuzz(
      eset = eset_i,
      centers = num_clusters_i,
      m = m_i
    )

    # Reorder clusters so those with similar trajectories are consecutive
    FCM_i <- .reorder_clusters(FCM_i)

    # Include input matrix (both modalities) and value of weighting exponent, m
    FCM_i[["input"]] <- exprs(eset_list[[tissue_i]])
    FCM_i[["call"]][["m"]] <- m_i

    cmeans_list[[tissue_i]] <- FCM_i
  }

  return(cmeans_list)
}



## Internal functions ----------------------------------------------------------

#' @title Reorder Fuzzy C-Means Clusters
#'
#' @description Reorder fuzzy c-means clusters so those with similar
#'   trajectories are consecutive.
#'
#' @param fclust an object of class \code{"fclust"}.
#'
#' @returns An object of class \code{"fclust"} with the clusters reordered and
#'   relabeled.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom stats cor as.dist hclust
#'
#' @noRd

.reorder_clusters <- function(fclust) {
  d <- as.dist(1 - cor(t(fclust[["centers"]])))
  hc <- hclust(d)
  neworder <- hc[["order"]]

  fclust[["centers"]] <- fclust[["centers"]][neworder, ]
  rownames(fclust[["centers"]]) <- seq_len(nrow(fclust[["centers"]]))

  fclust[["size"]] <- fclust[["size"]][neworder]

  fclust[["cluster"]][] <- match(fclust[["cluster"]], neworder)
  fclust[["membership"]] <- fclust[["membership"]][, neworder]
  colnames(fclust[["membership"]]) <- seq_len(ncol(fclust[["membership"]]))

  return(fclust)
}
