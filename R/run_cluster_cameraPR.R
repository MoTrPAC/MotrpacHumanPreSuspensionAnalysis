#' @title Analyze Fuzzy C-Means Clustering Results with CAMERA-PR
#'
#' @description Test for localization of molecular signatures to cluster
#'   centroids. Performs modified upper-tail signed rank tests with
#'   \code{TMSig::cameraPR.matrix} on the matrices of cluster membership
#'   probabilities obtained from fuzzy c-means (FCM) clustering. Identifies
#'   molecular signatures that closely follow the trajectories of each cluster's
#'   centroid.
#'
#' @param FCM fuzzy c-means (FCM) clustering results. Output of
#'   \code{\link{run_cmeans}}. This should be a named list of objects of class
#'   \code{fclust}.
#' @inheritParams run_cameraPR
#'
#' @returns An object of class \code{data.frame} with the following columns:
#'   \describe{
#'     \item{\code{tissue}}{factor; the tissue.}
#'     \item{\code{assay}}{factor; the omics assay.}
#'     \item{\code{cluster}}{factor; the cluster number. The clusters have the
#'     same meaning across omics assays measured in the same tissue, but they
#'     have no such relationship across tissues.}
#'     \item{\code{collection}}{factor; the broad molecular signature
#'     collection. See
#'     \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}} for details.}
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
#'     \item{\code{p_value}}{numeric; the upper-tail p-value.}
#'     \item{\code{adj_p_value}}{numeric; the BH-adjusted p-value. P-values are
#'     adjusted within each combination of tissue, assay, collection, and
#'     cluster.}
#'   }
#'
#' @seealso \code{\link{run_cmeans}}, \code{\link{run_cluster_ORA}},
#'   \code{\link[MotrpacHumanPreSuspensionAnalysis]{MOLECULAR_SIGNATURES}},
#'    \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}}
#'
#' @author Tyler Sagendorf
#'
#' @importFrom dplyr %>% bind_rows mutate select everything left_join rename across
#' @importFrom stats p.adjust
#'
#' @export run_cluster_cameraPR
#'
#' @examples
#' \dontrun{
#'   FCM <- run_cmeans()
#'
#'   # Run CAMERA-PR with all available molecular signatures
#'   cluster_res <- run_cluster_cameraPR(FCM = FCM)
#'   head(cluster_res)
#' }

run_cluster_cameraPR <- function(FCM,
                                 selected_omes = c("transcript-rna-seq",
                                                   "prot-pr", "prot-ph",
                                                   "metab"),
                                 selected_tissues = "all",
                                 database = names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
                                 path_to_gmt = NULL,
                                 min_size = 5L,
                                 overlap_cutoff = 0.7)
{
  on.exit(gc())

  # Allow users with older R versions to still use the package
  check_package_installation(pkg = "TMSig", fun = "run_cluster_cameraPR")

  # List of membership probability matrices by tissue and ome
  mem_list <- .prepare_cluster_mem(FCM = FCM,
                                   selected_omes = selected_omes,
                                   selected_tissues = selected_tissues)

  # Features that contributed to FCM results
  background_list <- lapply(mem_list, rownames)

  ## Prepare molecular signatures ----
  index <- .create_index(database = database,
                         path_to_gmt = path_to_gmt)

  index_list <- .prepare_sets(background_list = background_list,
                              index = index,
                              overlap_cutoff = overlap_cutoff,
                              min_size = min_size)

  res_list <- lapply(names(index_list), function(name_i) {
    TMSig::cameraPR.matrix(
      statistic = mem_list[[name_i]],
      index = index_list[[name_i]],
      use.ranks = TRUE, # modified signed rank test
      inter.gene.cor = 0.01,
      alternative = "greater", # upper-tail test
      adjust.globally = FALSE,
      min.size = min_size
    )
  })
  names(res_list) <- names(index_list)

  ## Reformat results ----
  out <- res_list %>%
    bind_rows(.id = "idcol") %>%
    dplyr::mutate(tissue = sub("\\..*", "", idcol),
           assay = sub(".*\\.", "", idcol),
           idcol = NULL,
           across(.cols = c(tissue, assay),
                  .fns = ~ factor(.x, levels = unique(.x)))) %>%
    dplyr::select(-Direction) %>% # all "Up"
    # Rename columns
    dplyr::rename(cluster = Contrast, set = GeneSet, set_size = NGenes,
           p_value = PValue, adj_p_value = FDR) %>%
    # Include collection, database, set_id, and set_short columns
    dplyr::left_join(MotrpacHumanPreSuspensionAnalysis::SET_TO_ID, by = "set") %>%
    droplevels.data.frame() %>%
    dplyr::mutate(set_size_DB = lengths(index)[set],
           size_ratio = round(set_size / set_size_DB, digits = 3L),
           across(.cols = everything(),
                  .fns = ~ structure(.x, names = NULL)),
           across(.cols = c(set_size, set_size_DB),
                  .fns = as.integer)) %>%
    # Convert set columns to factors to reduce the object size
    dplyr::mutate(across(.cols = c(set_id, set, set_short),
                  .fns = ~ factor(.x, levels = sort(unique(.x))))) %>%
    # Adjust p-values separately by tissue, ome, collection, and cluster.
    dplyr::mutate(.by = c(tissue, assay, collection, cluster),
           adj_p_value = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(tissue, assay, cluster, collection, database, p_value) %>%
    # Reorder columns
    dplyr::select(tissue, assay, cluster,
           collection, database, set_id, set, set_short,
           set_size, set_size_DB, size_ratio,
           p_value, adj_p_value) %>%
    # Remove columns with all NA values
    dplyr::select(where(function(x) !all(is.na(x))))

  return(out)
}



## Internal functions ----------------------------------------------------------

#' @title Extract cluster membership probability matrices from FCM results
#'
#' @inheritParams run_cluster_cameraPR
#'
#' @returns A named list of cluster membership probability matrices. Names are
#'   of the form "tissue.ome".
#'
#' @author Tyler Sagendorf
#'
#' @noRd

.prepare_cluster_mem <- function(FCM,
                                 selected_omes = c("transcript-rna-seq",
                                                   "prot-pr", "prot-ph",
                                                   "metab"),
                                 selected_tissues = "all")
{
  if (!is.vector(FCM, mode = "list") || is.null(names(FCM))) {
    stop("`FCM` must be the output of run_cmeans(). ",
         "See documentation for details.")
  }

  selected_omes <- match.arg(arg = tolower(selected_omes),
                             choices = c("transcript-rna-seq", "prot-pr",
                                         "prot-ph", "metab"),
                             several.ok = TRUE)

  selected_tissues <- match.arg(arg = tolower(selected_tissues),
                                choices = c("all", "adipose",
                                            "blood", "muscle"),
                                several.ok = TRUE)

  if (any(selected_tissues == "all")) {
    selected_tissues <- c("adipose", "blood", "muscle")
  }

  FCM <- FCM[names(FCM) %in% selected_tissues]

  if (!length(FCM)) {
    stop("names(FCM) do not match `selected_tissues`.")
  }

  mem_list <- lapply(FCM, function(fclust) {
    mem <- fclust[["membership"]]

    # Convert features to gene symbols, flanking sequences, or RefMet IDs and
    # split membership matrix rows by ome.
    mem_list <- .split_membership_by_ome(mem = mem,
                                         selected_omes = selected_omes)

    return(mem_list)
  })

  mem_list <- unlist(mem_list, recursive = FALSE)

  keep <- vapply(mem_list, function(mem_i) !is.null(mem_i), logical(1L))

  if (sum(keep) == 0L) {
    stop("Rownames of the membership probability matrices must ",
         "start with one of `selected_omes`.")
  }

  mem_list <- mem_list[keep]

  return(mem_list)
}


#' @title Split FCM Membership Probability Matrix by Ome and Convert Feature IDs
#'
#' @param mem cluster membership probability matrix. Row names must be of the
#'   form "assay feature_id".
#' @inheritParams .prepare_cluster_mem
#'
#' @returns A named list of membership probability matrices. Names are different
#'   omes.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom dplyr %>% filter select mutate case_when left_join arrange all_of desc
#' @importFrom tidyr pivot_longer pivot_wider nest unnest
#' @importFrom tibble rownames_to_column column_to_rownames deframe
#'
#' @noRd

.split_membership_by_ome <- function(mem, selected_omes) {
  mem_df <- mem %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_id") %>%
    dplyr::mutate(assay = sub("(^[^ ]+).*", "\\1", feature_id),
           feature_id = sub("[^ ]+ ", "", feature_id)) %>%
    dplyr::filter(assay %in% selected_omes)

  if (nrow(mem_df) == 0L)
    return(NULL)

  feature_conv <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(!grepl("^epi", assay)) %>%
    dplyr::select(feature_id, gene_symbol, flanking_sequence) %>%
    dplyr::mutate(across(.cols = everything(),
                  .fns = as.character)) %>%
    dplyr::mutate(
      new_id = case_when(
        !is.na(flanking_sequence) ~ flanking_sequence,
        !is.na(gene_symbol) ~ gene_symbol,
        TRUE ~ feature_id # if N/A or if features are RefMet names
      )
    ) %>%
    dplyr::mutate(new_id = strsplit(new_id, split = "\\|")) %>%
    tidyr::unnest(cols = new_id) %>%
    dplyr::select(feature_id, new_id)

  mem_list <- mem_df %>%
    dplyr::left_join(feature_conv, by = "feature_id") %>%
    # Some features are not in HUMAN_FEATURE_TO_GENE
    dplyr::mutate(new_id = ifelse(is.na(new_id), feature_id, new_id)) %>%
    tidyr::pivot_longer(cols = all_of(colnames(mem)),
                        names_to = "cluster",
                        values_to = "membership") %>%
    dplyr::mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
    # Select highest probability per new_id in each assay/cluster combination
    dplyr::arrange(assay, cluster, desc(membership), new_id) %>%
    dplyr::filter(.by = c(assay, cluster),
           !duplicated(new_id)) %>%
    tidyr::nest(.by = assay) %>%
    dplyr::mutate(data = lapply(data, function(data_i) {
      data_i %>%
        tidyr::pivot_wider(id_cols = new_id,
                           names_from = cluster,
                           values_from = membership) %>%
        tibble::column_to_rownames("new_id") %>%
        as.matrix()
    })) %>%
    tibble::deframe()

  return(mem_list)
}
