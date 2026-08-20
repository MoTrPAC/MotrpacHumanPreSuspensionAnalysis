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
#'   Naming any single metabolomics platform is the same as naming
#'   \code{"metab"}: the research platforms live in one stacked object per
#'   tissue, so all of them are loaded and the platform is read off the
#'   \code{platform} column. \code{"metab-t-clinical"} is the exception — it is
#'   clinical chemistry, is not in the stack, and is gated by
#'   \code{load_clinical}.
#'
#' @param single_matrix logical; if \code{TRUE}, returns a single combined
#'   \code{data.frame} across all selected tissues and assays. If \code{FALSE}
#'   (default), returns a nested list.
#'
#' @param verbose logical; toggle verbosity.
#'
#' @param load_clinical logical; whether to include the clinical chemistry omes
#'   (\code{clinical_ome_list()}: \code{"prot-clinical"} and
#'   \code{"metab-t-clinical"}). \code{FALSE} by default, so \code{"all"}
#'   returns the research omes and nothing changes for callers written before
#'   v2.0 split clinical chemistry out. Set \code{TRUE} to include them; they
#'   are dropped even when named unless it is set.
#'
#' @returns
#' If \code{single_matrix = FALSE}, a nested list of \code{data.frame} objects.
#' The top-level names correspond to tissues and the second-level names to
#' assays — the same nesting \code{\link{load_differential_analysis}} returns,
#' so the two tiers can be walked together. The research metabolomics platforms
#' arrive as a single \code{"metab"} element per tissue rather than one element
#' per platform.
#'
#' If \code{single_matrix = TRUE}, a single \code{data.frame} containing all
#' selected summary statistics, with missing columns filled as \code{NA}.
#'
#' Each table carries \code{tissue}, \code{assay}, \code{randomGroupCode},
#' \code{Timepoint}, \code{feature_id}, \code{Count}, \code{Mean} and \code{SD}.
#' The metabolomics tables carry \code{assay = "metab"} and one further column,
#' \code{platform}, naming the platform the row was measured on. That is how the
#' \code{*_DA} objects are labelled, so a join between the two tiers no longer
#' has to translate between two names for the same ome.
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
#' ## Load metabolomics only. One table per tissue, every platform in it.
#' metab_stats = load_summary_stats(selected_omes = "metab")
#' unique(metab_stats[["blood"]][["metab"]][["platform"]])
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
                              verbose = TRUE,
                              load_clinical = FALSE){

  selected_tissues <- match.arg(
    arg = selected_tissues,
    choices = c("all", "adipose", "blood", "muscle"),
    several.ok = TRUE
  )

  # Asking for one metabolomics platform loads them all, because there is one
  # object holding all of them: the research platforms are stacked into a single
  # {TISSUE}_METAB_SUM_STATS per tissue, keyed the way *_METAB_DA is keyed.
  #
  # metab-t-clinical is exempt. It is clinical chemistry, kept out of the stack
  # and published separately as BLOOD_METAB_T_CLINICAL_SUM_STATS; folding it into
  # "metab" would quietly return the stacked table instead of the one that was
  # asked for. It reaches selected_omes only when named, or through "all", and
  # the load_clinical gate below decides whether it survives either way.
  metab_platforms <- grepl("metab", selected_omes) &
    !selected_omes %in% clinical_ome_list()
  if (any(metab_platforms)) {
    selected_omes = c(selected_omes[!metab_platforms], "metab")
    if(verbose) message("By default, if any metab platform is loaded, all of them are loaded")
  }

  selected_omes <- match.arg(
    arg = selected_omes,
    choices = c(
      "all", "transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol",
      "epigen-atac-seq", "epigen-methylcap-seq", "metab",
      clinical_ome_list()
    ),
    several.ok = TRUE
  )

  if ("all" %in% selected_tissues) {
    selected_tissues <- c("adipose", "blood", "muscle")
  }

  if ("all" %in% selected_omes) {
    selected_omes <- c(
      "transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol",
      "epigen-atac-seq","epigen-methylcap-seq", "metab",
      clinical_ome_list()
    )
  }

  # Clinical chemistry is opt-in, the same way epigenomics is. Applied after
  # both expansions so it governs "all" and a named request alike.
  if (!load_clinical) {
    dropped <- base::intersect(selected_omes, clinical_ome_list())
    remaining <- base::setdiff(selected_omes, clinical_ome_list())
    # Asking only for what the gate removes leaves nothing to load, and an empty
    # selection surfaces further down as an error about a missing column. Say
    # what actually happened.
    if (length(dropped) && !length(remaining)) {
      stop("You've requested only clinical omes (",
           paste(dropped, collapse = ", "),
           ") but `load_clinical = FALSE`. Set `load_clinical = TRUE` to load ",
           "clinical chemistry.")
    }
    selected_omes <- remaining
    if (length(dropped) && verbose) {
      message("Clinical omes (", paste(dropped, collapse = ", "),
              ") are skipped; set `load_clinical = TRUE` to include them.")
    }
  }

  if("prot-ph" %in% selected_omes & verbose){
    message("Only features qualifying for diffential analysis are included. For proteomics and phosphoproteomics, this means some samples with missingness patterns that lead to paired n < 3 for any group are not included here.")
  }
  if(any(grepl("epigen", selected_omes)) & verbose){
    message("Epigenetics summary stats are trimmed to only show significant features due to file size limitations")
  }
  if("metab" %in% selected_omes & verbose){
    message("Metabolomics platforms are returned stacked in one table per tissue, with the
            platform in the `platform` column and `assay` reading \"metab\", matching the
            differential analysis. Please remember that the lowest CV metabolite is chosen and
            the relevant refmet name is used. If you're not able to find your desired
            metabolite, look through the METABOLOMICS_CVS object for the relevant
            refmet/feature name.")
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
