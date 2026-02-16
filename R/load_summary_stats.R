#' @title Load Summary Statistics for Normalized Expression Data
#'
#' @description
#' Loads group- and timepoint-level summary statistics for normalized expression
#' data across tissues, molecular assays, and analytical platforms. Summary
#' statistics consist of means and standard deviations computed within
#' randomization groups and timepoints. Limited to the initial acute bout, training
#' data is excluded here.
#'
#' These datasets are intended for descriptive and exploratory analyses.
#' Sample-level data are not distributed with this package and are available
#' upon request via \url{https://motrpac-data.org}.
#'
#' For metabolomics assays, summary statistics are computed after filtering
#' redundant metabolites; see the Methods section of the associated
#' documentation for details.
#'
#' @details
#' Summary statistics were filtered to only those that qualified for differential analysis.
#' This means for proteomics/phosphoproteomics, samples required a paired n>=3 to be included.
#' See the methods in the manuscript for more information.
#' In epigenetic assays, features were filtered for significant features (FDR<0.05) for file size purposes.
#'
#' @param selected_tissues character; tissues to include. One or more of
#'   \code{"adipose"}, \code{"blood"}, \code{"muscle"}, or \code{"all"}.
#'
#' @param selected_omes character; molecular assays to include. One or more of
#'   \code{"transcript-rna-seq"}, \code{"prot-pr"}, \code{"prot-ph"},
#'   \code{"prot-ol"}, \code{"epigen-atac-seq"},
#'   \code{"epigen-methylcap-seq"}, \code{"metab"}, or \code{"all"}.
#'   Selecting \code{"metab"} loads all metabolomics platforms.
#'
#' @param single_matrix logical; if \code{TRUE}, returns a single combined
#'   \code{data.frame} across all selected tissues and assays. If \code{FALSE}
#'   (default), returns a nested list.
#'
#' @param verbose logical; toggle verbosity.
#'
#' @returns
#' If \code{single_matrix = FALSE}, a nested list of \code{data.frame} objects.
#' The top-level names correspond to tissues, and the second-level names
#' correspond to assays or platforms.
#'
#' If \code{single_matrix = TRUE}, a single \code{data.frame} containing all
#' selected summary statistics, with missing columns filled as \code{NA}.
#'
#' @author Christopher Jin
#'
#' @export load_summary_stats
#'
#' @examples
#' \dontrun{
#' ## Load all summary statistics
#' sum_stats = load_summary_stats()
#'
#' ## Load metabolomics only
#' metab_stats = load_summary_stats(selected_omes = "metab")
#'
#' ## Load adipose transcriptomics as a single table
#' adipose_rna = load_summary_stats(
#'   selected_tissues = "adipose",
#'   selected_omes = "transcript-rna-seq",
#'   single_matrix = TRUE
#' )
#' }

load_summary_stats = function(selected_tissues = "all",
                              selected_omes = "all",
                              single_matrix = FALSE,
                              verbose = TRUE){

  selected_tissues <- match.arg(
    arg = selected_tissues,
    choices = c("all", "adipose", "blood", "muscle"),
    several.ok = TRUE
  )

  if(any(grepl("metab", selected_omes))){
    selected_omes = selected_omes[-grep("metab", selected_omes)]
    selected_omes = c(selected_omes, metab_only_list())
    if(verbose) message("By default, if any metab platform is loaded, all of them are loaded")
  }

  selected_omes <- match.arg(
    arg = selected_omes,
    choices = c(
      "all", "transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol",
      "epigen-atac-seq", "epigen-methylcap-seq", metab_only_list()
    ),
    several.ok = TRUE
  )

  if ("all" %in% selected_tissues) {
    selected_tissues <- c("adipose", "blood", "muscle")
  }

  if ("all" %in% selected_omes) {
    selected_omes <- c(
      "transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol",
      "epigen-atac-seq","epigen-methylcap-seq", metab_only_list()
    )
  }

  if("prot-ph" %in% selected_omes & verbose){
    message("Only features qualifying for diffential analysis are included. For proteomics and phosphoproteomics, this means some samples with missingness patterns that lead to paired n < 3 for any group are not included here.")
  }
  if(any(grepl("epigen", selected_omes)) & verbose){
    message("Epigenetics summary stats are trimmed to only show significant features due to file size limitations")
  }


  sum_stat_files <- data(package = "MotrpacHumanPreSuspensionAnalysis")
  sum_stat_files <- sum_stat_files[["results"]][, "Item"]
  sum_stat_files <- sum_stat_files[grepl("_SUM_STATS$", sum_stat_files)]

  tissues <- tolower(sub("\\_.*", "", sum_stat_files))

  omes <- sub("^[^_]+_(.*)_SUM_STATS$", "\\1", sum_stat_files)
  omes <- gsub("_", "-", tolower(omes))

  new_names <- structure(
    .Data = paste0(tissues, ".", omes),
    names = sum_stat_files
  )

  keep <- (tissues %in% selected_tissues) & (omes %in% selected_omes)
  new_names <- new_names[keep]

  # Load QC results into a list. Only works because of lazy loading
  out <- vector(mode = "list", length = length(new_names))
  names(out) <- as.character(new_names)

  for (i in seq_along(new_names)) {
    out[[i]] <- eval(parse(text = names(new_names[i])))
  }

  if (single_matrix) {
    out <- rbindlist(l = out, use.names = TRUE, fill = TRUE)
    setorderv(x = out, cols = "tissue", order = 1L)

  } else {

    tissues <- sub("\\..*$", "", names(out))
    names(out) <- sub(".*\\.", "", names(out))
    out <- split(do.call(list, out), tissues)
  }

  return(out)

}
