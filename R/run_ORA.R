#' @title Over-Representation Analysis (ORA)
#'
#' @description Fast over-representation analysis.
#'
#' @param input character; vector of "interesting" features. Most likely genes,
#'   but may be RefMet metabolite IDs or singly-phosphorylated peptides.
#' @param background character; vector of all features that were analyzed. Used
#'   to filter the molecular signatures.
#' @inheritParams run_cameraPR
#'
#' @returns An object of class \code{data.frame} with the following columns:
#'   \describe{
#'     \item{\code{collection}}{factor; the broad molecular signature
#'     collection. Only included when \code{path_to_gmt} is \code{NULL}. See
#'     \code{\link{SET_TO_ID}} for details.}
#'     \item{\code{database}}{factor; the molecular signature database. Only
#'     included when \code{path_to_gmt} is \code{NULL}. See
#'     \code{\link{SET_TO_ID}} for details.}
#'     \item{\code{set_id}}{character; a unique ID for the molecular
#'     signature. Only included when \code{path_to_gmt} is \code{NULL}. See
#'     \code{\link{SET_TO_ID}} for details.}
#'     \item{\code{set}}{character; the molecular signature being tested.
#'     For global proteomics and transcriptomics, these are gene sets. For
#'     phosphoproteomics, these are kinase sets.}
#'     \item{\code{set_short}}{character; a shortened version of
#'     \code{set}. Only included when \code{path_to_gmt} is \code{NULL}. See
#'     \code{\link{SET_TO_ID}} for details.}
#'     \item{\code{set_size}}{integer; the number of molecules in the set that
#'     were present in the \code{background} vector.}
#'     \item{\code{set_size_DB}}{integer; the number of molecules in the set,
#'     as defined in \code{\link{MOLECULAR_SIGNATURES}}.}
#'     \item{\code{size_ratio}}{numeric; the ratio of \code{set_size} to
#'     \code{set_size_DB}, rounded to the nearest thousandth. A measure of
#'     confidence that the gene set being tested is correctly described by the
#'     entry in the \code{set} column. While smaller values do not
#'     necessarily indicate that the results are unreliable, terms from the gene
#'     set databases should be treated with caution.}
#'     \item{\code{set_size_in_input}}{integer; the number of molecules in the
#'     set that were present in the \code{input} vector.}
#'     \item{\code{input_size}}{integer; the size of the \code{input} vector.}
#'     \item{\code{background_size}}{integer; the size of the \code{background}
#'     vector.}
#'     \item{\code{p_value}}{numeric; the two-sided p-value.}
#'     \item{\code{adj_p_value}}{numeric; the BH-adjusted p-value. P-values
#'     are adjusted within each combination of tissue, assay, contrast, and
#'     collection.}
#'   }
#'
#' @seealso \code{\link{run_cluster_ORA}}, \code{\link{run_cameraPR}}
#'
#' @importFrom dplyr %>% mutate across arrange left_join select any_of everything
#' @importFrom stats phyper p.adjust
#'
#' @export run_ORA
#'
#' @examples
#' # Use genes from all gene sets as the background
#' bg <- MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES
#' bg[c("PSP", "REFMET")] <- NULL
#' bg <- unique(unlist(bg))
#'
#' # Use all genes from the MitoCarta OXPHOS term as the interesting subset
#' input <- MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES$MITOCARTA[["MITOCARTA_OXPHOS"]]
#'
#' res <- run_ORA(input = input, background = bg)
#' head(res)

run_ORA <- function(input,
                    background,
                    database = names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
                    path_to_gmt = NULL,
                    min_size = 5L,
                    overlap_cutoff = 0.7) {
  on.exit(gc())

  # Allow users with older R versions to still use the package
  check_package_installation(pkg = "TMSig", fun = "run_ORA")

  overlap_cutoff <- max(0, min(overlap_cutoff, 1, na.rm = TRUE))

  if (!is.vector(input, mode = "character") || !length(input)) {
    stop("`input` must be a character vector of features.")
  }

  if (!is.vector(background, mode = "character") || !length(background)) {
    stop("`background` must be a character vector of features.")
  }

  input <- unique(input[!is.na(input)])
  background <- unique(background[!is.na(background)])

  if (any(!input %in% background)) {
    stop("`input` must be a subset of `background`.")
  }

  ## Prepare molecular signatures ----
  index <- .create_index(database = database,
                         path_to_gmt = path_to_gmt)

  set_size_DB <- lengths(index)

  if (all(grepl("^PTMSIGDB", names(index)))) {
    # This code is from TMSig::filterSets. It was repurposed to work with
    # PTMsigDB, where the sites in each set end with ";u" or ";d", but the
    # background vector of IDs do not.
    set_dt <- data.table(
      sets = rep(names(index), lengths(index)),
      elements = unlist(index),
      stringsAsFactors = FALSE
    )

    set_dt[, elements2 := sub(";.*$", "", elements)]
    set_dt <- subset(set_dt, subset = elements2 %in%  background)

    index <- split(x = set_dt[["elements2"]], f = set_dt[["sets"]])
    set_sizes <- lengths(index)
    keep_sizes <- (set_sizes >= min_size) & (set_sizes < length(background))
    index <- index[keep_sizes]
  } else {
    # Restrict sets to background, filter by size
    index <- TMSig::filterSets(
      x = index,
      background = background,
      min_size = min_size,
      max_size = length(background) - 1L
    )
  }

  # overlap_cutoff does not affect RefMet, PhosphoSitePlus, or PTMSigDB
  # signatures
  if (!any(grepl("^REFMET|^PSP|^PTMSIGDB", names(index)))) {
    keep <- lengths(index) / set_size_DB[names(index)] >= overlap_cutoff

    if (sum(keep) == 0L) {
      stop("No gene sets pass `overlap_cutoff`. ",
           "The `background` may be too small.")
    }

    index <- index[keep]
  }

  # Restrict sets to input vector of features.
  index_filt <- .fast_list_intersect(x = index, y = input)
  # New chnge as of Oct 8th 2025: empty sets are NOT kept.
  #so we dont test any pathways without a matching value.
  index_filt = index_filt[lengths(index_filt) > 0]


  out <- data.frame(set = names(index),
                    set_size = lengths(index)) %>%
    dplyr::mutate(set_size_DB = set_size_DB[set],
           size_ratio = round(set_size / set_size_DB, digits = 3L),
           set_size_in_input = lengths(index_filt)[set],
           input_size = length(input),
           background_size = length(background),
           p_value = phyper(q = set_size_in_input - 1L,
                            m = set_size,
                            n = background_size - set_size,
                            k = input_size,
                            lower.tail = FALSE),
           across(.cols = everything(),
                  .fns = ~ structure(.x, names = NULL))) %>%
    dplyr::arrange(p_value) %>%
    dplyr::left_join(MotrpacHumanPreSuspensionAnalysis::SET_TO_ID, by = "set") %>%
    dplyr::mutate(.by = collection,
           adj_p_value = p.adjust(p_value, method = "BH")) %>%
    # Reorder columns
    dplyr::select(any_of(colnames(MotrpacHumanPreSuspensionAnalysis::SET_TO_ID)), everything()) %>%
    dplyr::filter(!is.na(p_value))

  return(out)
}
