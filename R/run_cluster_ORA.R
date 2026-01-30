#' @title Analyze Fuzzy C-means Clustering Results with ORA
#'
#' @description Over-representation analysis (ORA) applied to results of fuzzy
#'   c-means (FCM) clustering.
#'
#' @inheritParams run_cluster_cameraPR
#' @param min_prob numeric; the minimum cluster membership probability required
#'   for features to belong to a cluster. All features will be used for the
#'   background for ORA, regardless of membership probability.
#'
#' @returns An object of class \code{data.frame} with the following columns:
#'   \describe{
#'     \item{\code{tissue}}{factor; the tissue.}
#'     \item{\code{assay}}{factor; the omics assay.}
#'     \item{\code{cluster}}{factor; the cluster number. The clusters have the
#'     same meaning across omics assays measured in the same tissue, but they
#'     have no such relationship across tissues.}
#'     \item{\code{collection}}{factor; the broad molecular signature
#'     collection. See \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}}
#'     for details.}
#'     \item{\code{database}}{factor; the molecular signature database. See
#'     \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}} for details.}
#'     \item{\code{set_id}}{character; a unique ID for the molecular
#'     signature. See \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}}
#'     for details.}
#'     \item{\code{set}}{character; the molecular signature being tested.
#'     For global proteomics and transcriptomics, these are gene sets. For
#'     phosphoproteomics, these are kinase sets.}
#'     \item{\code{set_short}}{character; a shortened version of
#'     \code{set}. See \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}}
#'     for details.}
#'     \item{\code{set_size}}{integer; the number of molecules in the set that
#'     were present in the DA results for that specific tissue/assay
#'     combination.}
#'     \item{\code{set_size_DB}}{integer; the number of molecules in the set,
#'     as defined in the GMT file.}
#'     \item{\code{size_ratio}}{numeric; the ratio of \code{set_size} to
#'     \code{set_size_DB}, rounded to the nearest thousandth. A measure of
#'     confidence that the gene set being tested is correctly described by the
#'     entry in the \code{set} column. While smaller values do not
#'     necessarily indicate that the results are unreliable, terms from the gene
#'     set databases should be treated with caution.}
#'     \item{\code{set_size_in_cluster}}{integer; the number of elements in the
#'     set that localized to the cluster.}
#'     \item{\code{cluster_size}}{integer; number of features in each cluster
#'     with membership probabilities of at least \code{min_prob}.}
#'     \item{\code{background_size}}{integer; the number of elements in the
#'     background. Includes all features that appeared in clustering results.}
#'     \item{\code{p_value}}{numeric; the upper-tail p-value.}
#'     \item{\code{adj_p_value}}{numeric; the BH-adjusted p-value. P-values are
#'     adjusted within each combination of tissue, assay, collection, and
#'     cluster.}
#'   }
#'
#' @seealso \code{\link{run_cmeans}}, \code{\link{run_cluster_cameraPR}},
#'   \code{\link[MotrpacHumanPreSuspensionAnalysis]{MOLECULAR_SIGNATURES}},
#'   \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}}
#'
#' @author Tyler Sagendorf
#'
#' @importFrom dplyr %>% bind_rows mutate select everything left_join across
#' @importFrom stats p.adjust
#'
#' @export run_cluster_ORA
#'
#' @examples
#' \dontrun{
#'   FCM <- run_cmeans()
#'
#'   # Run ORA with all available molecular signatures
#'   cluster_res <- run_cluster_ORA(FCM = FCM)
#'   head(cluster_res)
#' }

run_cluster_ORA <- function(FCM,
                            selected_omes = c("transcript-rna-seq",
                                              "prot-pr", "prot-ph",
                                              "metab"),
                            selected_tissues = "all",
                            database = names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
                            path_to_gmt = NULL,
                            min_size = 5L,
                            overlap_cutoff = 0.7,
                            min_prob = 0.3)
{
  on.exit(gc())

  # Allow users with older R versions to still use the package
  check_package_installation(pkg = "TMSig", fun = "run_cluster_ORA")

  min_prob <- max(0, min(1, min_prob, na.rm = TRUE))

  # List of membership probability matrices by tissue and ome
  mem_list <- .prepare_cluster_mem(FCM = FCM,
                                   selected_omes = selected_omes,
                                   selected_tissues = selected_tissues)

  # Features in any cluster, regardless of membership probability
  background_list <- lapply(mem_list, rownames)

  # Prepare list of clusters. Features are filtered based on membership
  # probability first
  cluster_list <- lapply(mem_list, function(mem_i) {
    cluster_levels <- seq_len(ncol(mem_i))

    keep <- apply(mem_i >= min_prob, 1, any)
    mem_i <- mem_i[keep, , drop = FALSE]
    cluster_id <- apply(mem_i, 1, which.max)

    # Convert to factor to keep empty clusters
    cluster_id <- factor(cluster_id, levels = cluster_levels)

    split(x = rownames(mem_i), f = cluster_id)
  })

  ## Prepare molecular signatures ----
  index <- .create_index(database = database,
                         path_to_gmt = path_to_gmt)

  index_list <- .prepare_sets(background_list = background_list,
                              index = index,
                              overlap_cutoff = overlap_cutoff,
                              min_size = min_size)

  ## Over-representation analysis ----
  res_list <- lapply(names(index_list), function(name_i) {
    .cluster_ORA(
      cluster_list = cluster_list[[name_i]],
      background = background_list[[name_i]],
      index = index_list[[name_i]],
      min_size = min_size
    )
  })
  names(res_list) <- names(index_list)

  ## Reformat results ----
  out <- res_list %>%
    bind_rows(.id = "idcol") %>%
    mutate(tissue = sub("\\..*", "", idcol),
           assay = sub(".*\\.", "", idcol),
           idcol = NULL,
           across(.cols = c(tissue, assay),
                  .fns = ~ factor(.x, levels = unique(.x)))) %>%
    # Include collection, database, set_id, and set_short columns
    left_join(MotrpacHumanPreSuspensionAnalysis::SET_TO_ID, by = "set") %>%
    droplevels.data.frame() %>%
    mutate(set_size_DB = lengths(index)[set],
           size_ratio = round(set_size / set_size_DB, digits = 3L),
           across(.cols = everything(),
                  .fns = ~ structure(.x, names = NULL)),
           across(.cols = c(set_size, set_size_DB),
                  .fns = as.integer)) %>%
    # Convert set columns to factors to reduce the object size
    mutate(across(.cols = c(set_id, set, set_short),
                  .fns = ~ factor(.x, levels = sort(unique(.x))))) %>%
    # Adjust p-values separately by tissue, ome, collection, and cluster.
    mutate(.by = c(tissue, assay, collection, cluster),
           adj_p_value = p.adjust(p_value, method = "BH")) %>%
    arrange(tissue, assay, cluster, collection, database, p_value) %>%
    # Reorder columns
    select(tissue, assay, cluster,
           collection, database, set_id, set, set_short,
           set_size, set_size_DB, size_ratio,
           set_size_in_cluster, cluster_size, background_size,
           p_value, adj_p_value) %>%
    # Remove columns with all NA values
    select(where(function(x) !all(is.na(x))))

  return(out)
}



## Internal functions ----------------------------------------------------------

#' @title Over-representation analysis of hard clustering results
#'
#' @description Over-representation analysis of hard clustering results.
#'
#' @param cluster_list a named list of elements in each cluster. Elements must
#'   be of type \code{character}.
#' @inheritParams run_cameraPR
#'
#' @returns An object of class \code{data.frame} containing the ORA results.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom data.table data.table := rbindlist setorderv setcolorder setDF
#' @importFrom stats phyper p.adjust
#'
#' @noRd

.cluster_ORA <- function(cluster_list,
                         background,
                         index,
                         min_size = 5L)
{
  background <- unique(background)
  background <- background[!is.na(background)]
  size_background <- length(background)

  index <- TMSig::filterSets(
    x = index,
    background = background,
    min_size = min_size,
    max_size = size_background
  )
  set_sizes <- lengths(index)

  dt_list <- lapply(cluster_list, function(cluster_i) {
    cluster_i_elements <- unlist(cluster_i)
    size_cluster_i <- length(cluster_i_elements)

    # Restrict index to cluster elements. Keep empty sets
    index_ci <- .fast_list_intersect(x = index, y = cluster_i_elements)
    set_sizes_ci <- lengths(index_ci)

    dt_i <- data.table(set = names(set_sizes_ci),
                       set_size_in_cluster = set_sizes_ci,
                       stringsAsFactors = FALSE)

    return(dt_i)
  })

  dt <- rbindlist(dt_list, idcol = "cluster")

  dt[, `:=`(set_size = set_sizes[set],
            cluster_size = lengths(cluster_list)[cluster],
            background_size = size_background)]

  dt[, p_value := phyper(q = set_size_in_cluster - 1L,   # successes in cluster
                         m = set_size,                   # total successes
                         n = background_size - set_size, # total failures
                         k = cluster_size,               # sample_size
                         lower.tail = FALSE)]

  # Adjust p-values across clusters and sets
  dt[, adj_p_value := p.adjust(p_value, method = "BH")]
  dt[, cluster := factor(cluster, levels = names(cluster_list))]

  setorderv(dt, cols = c("cluster", "p_value", "set_size", "set"),
            order = c(1, 1, -1, 1))
  setcolorder(dt, neworder = "p_value", before = "adj_p_value")
  setDF(dt)

  return(dt)
}


#' @title Fast intersection of a list with a character vector
#'
#' @description Fast intersection of a list with a character vector. Much faster
#'   than \code{lapply(x, function(xi) intersect(xi, y))} if \code{lengths(x)}
#'   are small.
#'
#' @param x a named list of character vectors.
#' @param y a character vector. Missing values and duplicates will be removed.
#'
#' @returns A named list with the same size as x, unless \code{names(x)} are
#'   duplicated.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom data.table data.table := setorderv
#'
#' @noRd

# Note: does not check that x is a valid named list of character vectors or that
# y is a character vector.
.fast_list_intersect <- function(x, y) {
  # Convert list to data.table for fast filtering
  dt <- data.table(sets = rep(names(x), lengths(x)),
                   elements = unlist(x),
                   stringsAsFactors = FALSE)

  # Convert to factor to preserve order and keep empty sets when splitting
  dt[, sets := factor(sets, levels = unique(names(x)))]

  setorderv(dt, cols = c("sets", "elements"), order = rep(1L, 2L))

  y <- unique(y)
  y <- y[!is.na(y)]
  dt <- subset(dt, subset = elements %in% y)
  dt <- unique(dt)

  out <- split(x = dt[["elements"]], f = dt[["sets"]])

  return(out)
}

