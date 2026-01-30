#' @title Wrapper for CAMERA-PR
#'
#' @description A wrapper for \code{\link[TMSig]{cameraPR.matrix}} that performs
#'   pre-ranked Correlation Adjusted MEan RAnk molecular signature analysis of
#'   differential analysis results.
#'
#' @param DA_list list; a named list of \code{data.frame} objects, each
#'   containing differential analysis results for a specific tissue/assay
#'   combination. The list must be nested, with tissues at the top level and
#'   omes within tissues. If \code{NULL} (default), the differential analysis
#'   results will be generated with
#'   \code{\link[MotrpacHumanPreSuspensionAnalysis]{load_differential_analysis}}.
#'   Unless wishing to analyze DA results that are not in
#'   MotrpacHumanPreSuspensionAnalysis, this should remain \code{NULL}.
#' @param selected_omes character; one or more character strings selected from
#'   the following options: \code{"transcript-rna-seq"}, \code{"prot-pr"},
#'   \code{"prot-ph"}, and \code{"metab"} (all metabolomics platforms). Passed
#'   to \code{\link[MotrpacHumanPreSuspensionAnalysis]{load_differential_analysis}}.
#' @param selected_tissues character; passed to
#'   \code{\link[MotrpacHumanPreSuspensionAnalysis]{load_differential_analysis}}.
#'   One or more of the following: \code{"all"}, \code{"muscle"},
#'   \code{"adipose"}, or \code{"blood"}.
#' @param database character; one or more names specifying the database(s) to
#'   test. Options are (case insensitive) \code{"BIOCARTA"},
#'   \code{"KEGG_MEDICUS"}, \code{"PID"}, \code{"REACTOME"}, \code{"WP"}
#'   (WikiPathways database), \code{"GOBP"}, \code{"GOCC"}, \code{"GOMF"},
#'   \code{"MITOCARTA"} (MitoCarta3.0 database), \code{"PSP"} (PhosphoSitePlus
#'   kinases; only valid when \code{selected_omes} contains \code{"prot-ph"}),
#'   or \code{"REFMET"} (RefMet chemical subclasses; only valid when
#'   \code{selected_omes} contains \code{"metab"}). See
#'   \code{\link[MotrpacHumanPreSuspensionAnalysis]{MOLECULAR_SIGNATURES}} for
#'   details.
#' @param path_to_gmt character; (optional) path to one or more GMT files.
#'   Passed to \code{TMSig::readGMT}. If provided, \code{database} is ignored.
#' @param min_size integer; the minimum set size for testing.
#' @param overlap_cutoff numeric; the minimum proportion of genes in each set
#'   that must appear in a given dataset. Used to pre-filter sets. Does not
#'   affect \code{"metab"} or \code{"prot-ph"} results. This will always be 0.1
#'   for \code{"prot-ol"} results.
#'
#' @returns An object of class \code{data.frame} with the following columns:
#'   \describe{
#'     \item{\code{tissue}}{factor; the tissue.}
#'
#'     \item{\code{assay}}{factor; the omics assay.}
#'
#'     \item{\code{contrast_type}}{factor; the type of contrast.}
#'
#'     \item{\code{contrast}}{factor; the contrast of interest.}
#'
#'     \item{\code{contrast_short}}{factor; shortened contrasts.}
#'
#'     \item{\code{collection}}{factor; the broad molecular signature
#'     collection. Only included when \code{path_to_gmt} is \code{NULL}. See
#'     \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}} for details.}
#'
#'     \item{\code{database}}{factor; the molecular signature database. Only
#'     included when \code{path_to_gmt} is \code{NULL}. See
#'     \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}} for details.}
#'
#'     \item{\code{set_id}}{character; a unique ID for the molecular
#'     signature. See \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}}
#'     for details.}
#'
#'     \item{\code{set}}{character; the molecular signature being tested.
#'     For global proteomics and transcriptomics, these are gene sets. For
#'     phosphoproteomics, these are kinase sets.}
#'
#'     \item{\code{set_short}}{character; a shortened version of
#'     \code{set}. Only included when \code{path_to_gmt} is \code{NULL}. See
#'     \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}} for details.}
#'
#'     \item{\code{set_size}}{integer; the number of molecules in the set that
#'     were present in the DA results for that specific tissue/assay
#'     combination.}
#'
#'     \item{\code{set_size_DB}}{integer; the number of molecules in the set,
#'     as defined in the GMT file.}
#'
#'     \item{\code{size_ratio}}{numeric; the ratio of \code{set_size} to
#'     \code{set_size_DB}, rounded to the nearest thousandth. A measure of
#'     confidence that the gene set being tested is correctly described by the
#'     entry in the \code{set} column. While smaller values do not
#'     necessarily indicate that the results are unreliable, terms from the gene
#'     set databases should be treated with caution.}
#'
#'     \item{\code{direction}}{factor; the direction of change: \code{"Up"}
#'     or \code{"Down"}.}
#'
#'     \item{\code{t}}{numeric; the two-sample t-statistic.}
#'
#'     \item{\code{df}}{integer; the available degrees of freedom.}
#'
#'     \item{\code{z.std}}{numeric; standard Normal equivalents of \code{t}.}
#'
#'     \item{\code{p_value}}{numeric; the two-sided p-value.}
#'
#'     \item{\code{adj_p_value}}{numeric; the BH-adjusted p-value. P-values
#'     are adjusted within each combination of tissue, assay, contrast, and
#'     collection.}
#'   }
#'
#' @seealso \code{\link[MotrpacHumanPreSuspensionAnalysis]{MOLECULAR_SIGNATURES}},
#'   \code{\link[MotrpacHumanPreSuspensionAnalysis]{SET_TO_ID}}
#'
#' @author Tyler Sagendorf
#'
#' @import MotrpacHumanPreSuspensionAnalysis
#' @importFrom dplyr %>% any_of bind_rows everything left_join mutate select
#'   rename across
#' @importFrom stats p.adjust
#' @importFrom utils modifyList
#'
#' @export run_cameraPR
#'
#' @examples
#' \dontrun{
#' # Test Reactome and Gene Ontology Biological Processes
#' # databases on muscle global proteomics.
#' x1 <- run_cameraPR(selected_omes = "prot-pr",
#'                    selected_tissues = "muscle",
#'                    database = c("GOBP", "REACTOME"))
#'
#' # Test all omes using all molecular signature databases
#' x2 <- run_cameraPR()
#' }

run_cameraPR <- function(DA_list = NULL,
                         selected_omes = c("transcript-rna-seq",
                                           "prot-pr",
                                           "prot-ph",
                                           "prot-ol",
                                           "metab"),
                         selected_tissues = "all",
                         database = setdiff(names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
                                            "PTMSIGDB"),
                         path_to_gmt = NULL,
                         min_size = 5L,
                         overlap_cutoff = 0.7)
{
  on.exit(gc())

  # Allow users with older R versions to still use the package
  check_package_installation(pkg = "TMSig", fun = "run_cameraPR")

  ## Prepare differential analysis results and molecular signatures ----
  ls <- .prepare_DA_results_and_sets(
    DA_list = DA_list,
    selected_omes = selected_omes,
    selected_tissues = selected_tissues,
    database = database,
    path_to_gmt = path_to_gmt,
    min_size = min_size,
    overlap_cutoff = overlap_cutoff
  )

  # Dump contents of ls into function namespace
  for (name_i in names(ls))
    assign(x = name_i, value = ls[[name_i]])

  ## CAMERA-PR ----
  res_list <- lapply(names(index_list), function(name_i) {
    TMSig::cameraPR.matrix(
      statistic = DA_list[[name_i]], # matrix of z-statistics
      index = index_list[[name_i]], # named list of signatures
      use.ranks = FALSE,
      inter.gene.cor = 0.01,
      sort = TRUE,
      alternative = "two.sided",
      adjust.globally = FALSE,
      min.size = min_size
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
                  .fns = ~ factor(.x, levels = sort(unique(.x))))) %>%
    # Rename columns
    rename(contrast = Contrast, set = GeneSet, set_size = NGenes,
           direction = Direction, t = TwoSampleT, z.std = ZScore,
           p_value = PValue, adj_p_value = FDR) %>%
    # Include contrast_type and contrast_short columns
    dplyr::left_join(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER, by = "contrast") %>%
    dplyr::mutate(
      contrast = factor(
        x = contrast,
        levels = intersect(
          levels(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER$contrast),
          levels(contrast)
        )),
      contrast_short = factor(
        x = contrast_short,
        levels = intersect(
          levels(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER$contrast_short),
          levels(contrast_short)
        ))
    ) %>%
    # Include collection, database, set_id, and set_short columns
    left_join(MotrpacHumanPreSuspensionAnalysis::SET_TO_ID, by = "set") %>%
    mutate(set_size = as.integer(set_size),
           set_size_DB = lengths(index)[set],
           size_ratio = round(set_size / set_size_DB,
                              digits = 3L),
           df = as.integer(df),
           direction = factor(direction, levels = c("Up", "Down")),
           across(.cols = everything(),
                  .fns = ~ structure(.x, names = NULL))) %>%
    # Convert set columns to factors to reduce the object size
    mutate(across(.cols = c(set_id, set, set_short),
                  .fns = ~ factor(.x, levels = sort(unique(.x))))) %>%
    droplevels.data.frame() %>%
    # Adjust p-values separately by tissue, ome, contrast, and collection
    mutate(.by = c(tissue, assay, contrast, collection),
           adj_p_value = p.adjust(p_value, method = "BH")) %>%
    arrange(contrast_type, tissue, assay, contrast,
            collection, database, p_value) %>%
    # Reorder columns
    select(tissue, assay, contrast_type, contrast, contrast_short,
           collection, database, set_id, set, set_short,
           set_size, set_size_DB, size_ratio, direction,
           t, df, z.std, p_value, adj_p_value)

  return(out)
}



## Helper functions ------------------------------------------------------------

#' @title Convert list of DA results tables to a list of z-statistic matrices
#'
#' @description Convert list of DA results to a list of z-statistic matrices
#'   with genes, RefMet metabolite IDs, or phosphosites as rows and contrasts as
#'   columns. Duplicates are resolved by selecting the most extreme z-score
#'   (positive or negative).
#'
#' @inheritParams run_cameraPR
#' @param convert_features logical; whether features should be converted to gene
#'   symbols, flanking sequences, or RefMet metabolite names. \code{TRUE} for
#'   enrichment analysis, \code{FALSE} for fuzzy c-means clustering.
#'
#' @returns a list of z-statistic matrices with features or gene symbols, RefMet
#'   metabolite IDs, or flanking sequences as rows and contrasts as columns.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom dplyr %>% bind_rows mutate left_join filter select arrange across
#'   where
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider unnest
#' @importFrom MotrpacHumanPreSuspensionAnalysis load_differential_analysis
#'
#' @noRd

.prepare_DA_results <- function(DA_list = NULL,
                                selected_omes = c("transcript-rna-seq",
                                                  "prot-pr",
                                                  "prot-ph",
                                                  "prot-ol",
                                                  "metab"),
                                selected_tissues = "all",
                                convert_features = TRUE,
                                .contrast_type = NULL)
{
  selected_omes <- match.arg(selected_omes,
                             choices = c("transcript-rna-seq", "prot-pr",
                                         "prot-ph", "metab", "prot-ol"),
                             several.ok = TRUE)

  selected_tissues <- match.arg(selected_tissues,
                                choices = c("all", "adipose",
                                            "blood", "muscle"),
                                several.ok = TRUE)

  if ("all" %in% selected_tissues) {
    selected_tissues <- c("adipose", "blood", "muscle")
  }

  if (is.null(DA_list)) {
    DA_list <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
      selected_omes = selected_omes,
      selected_tissues = selected_tissues,
      combine_with_featgene = FALSE
    )
  }

  # Combine metabolomics platforms
  # DA_list <- lapply(DA_list, function(tissue_i) {
  #   is_metab <- grep("metab", names(tissue_i))
  #
  #   if (length(is_metab)) {
  #     metab_df <- bind_rows(tissue_i[is_metab], .id = "platform")
  #
  #     tissue_i[is_metab] <- NULL
  #     tissue_i["metab"] <- list(metab_df)
  #   }
  #
  #   return(tissue_i)
  # })

  # Unnest list and collapse tissue and ome with a "."
  DA_list <- unlist(DA_list, recursive = FALSE)

  keep_tissues <- sub("\\..*", "", names(DA_list)) %in% selected_tissues
  keep_omes <- sub(".*\\.", "", names(DA_list)) %in% selected_omes

  DA_list <- DA_list[keep_tissues & keep_omes]
  DA_names <- names(DA_list)

  if (is.null(DA_names))
    stop("`DA_list` must be a named list of differential analysis results ",
         "tables with names of the form 'tissue.assay'.")

  # Used by run_cmeans
  if (!is.null(.contrast_type)) {
    DA_list <- lapply(DA_list, function(xi) {
      xi %>%
        select(-any_of(c("contrast_type", "contrast_short"))) %>%
        mutate(
          contrast = factor(contrast,
                            levels = levels(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER$contrast))
        ) %>%
        left_join(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER,
                  by = "contrast") %>%
        filter(contrast_type %in% .contrast_type) %>%
        arrange(contrast) %>%
        droplevels.data.frame()
    })
  }

  # Convert list of data.frames to a list of matrices with features (genes,
  # phosphosites, or metabolites/lipids) as rows and contrasts as columns.
  DA_list <- lapply(DA_names, function(name_i) {
    ome_i <- sub(".*\\.", "", name_i)

    x <- DA_list[[name_i]] %>%
      mutate(across(.cols = where(is.factor),
                    .fns = as.character))

    if (convert_features) {
      if (ome_i == "prot-ph") {
        flanking <-
          MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE[, c("feature_id", "flanking_sequence")] %>%
          mutate(across(.cols = everything(),
                        .fns = as.character))

        # Use flanking sequence as ID
        x <- x %>%
          left_join(flanking, by = "feature_id") %>%
          mutate(flanking_sequence = strsplit(flanking_sequence,
                                              split = "\\|")) %>%
          # Convert to single-sequence data
          tidyr::unnest(cols = flanking_sequence) %>%
          mutate(new_id = flanking_sequence)
      } else if (ome_i == "metab") {
        x <- x %>%
          mutate(feature_id = as.character(feature_id),
                 new_id = feature_id) %>%
          select(contrast, new_id, z.std)
      } else {
        # Used to convert proteins and transcripts to genes
        feature_to_symbol <-
          MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
          select(feature_id, gene_symbol) %>%
          mutate(across(.cols = where(is.factor),
                        .fns = as.character)) %>%
          distinct()

        ## "prot-pr", "transcript-rna-seq", and default behavior for
        ## user-supplied DA results that are not "prot-ph" or "metab".
        required_cols <- c("feature_id", "contrast", "z.std")

        if (any(!required_cols %in% colnames(x)))
          stop("All `DA_list` tables must contain the following columns: ",
               paste(dQuote(required_cols), collapse = ", "), ".")

        # Include gene_symbol column
        x <- left_join(x, feature_to_symbol,
                       by = "feature_id") %>%
          mutate(new_id = gene_symbol)

        if (all(is.na(x$gene_symbol)))
          stop("No features in the ", name_i,
               " DA results map to gene symbols.")
      }

      # For each gene (or phosphosite or metabolite/lipid) and contrast pairing,
      # select the most extreme z-score (positive or negative)
      x <- x %>%
        # If the new_id is missing, use the feature ID to avoid unnecessarily
        # collapsing or removing rows.
        mutate(new_id = ifelse(!is.na(new_id),
                               new_id,
                               feature_id)) %>%
        arrange(contrast, desc(abs(z.std)), z.std) %>%
        filter(.by = contrast,
               !duplicated(new_id))
    } else {
      # Do not convert features
      x <- mutate(x, new_id = feature_id)
    }

    # Convert to a matrix with contrasts as columns and features as rows. Values
    # of the matrix are z-scores.
    x_wide <- x %>%
      select(contrast, new_id, z.std) %>%
      tidyr::pivot_wider(id_cols = new_id,
                         names_from = contrast,
                         values_from = z.std) %>%
      tibble::column_to_rownames("new_id") %>%
      as.matrix()

    return(x_wide)
  })

  names(DA_list) <- DA_names

  return(DA_list)
}



#' @title Create a list of molecular signatures
#'
#' @inheritParams run_cameraPR
#'
#' @author Tyler Sagendorf
#'
#' @noRd

.create_index <- function(
    database = names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
    path_to_gmt = NULL
) {
  choices <- names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES)

  if (is.null(path_to_gmt)) {
    choices <- names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES)

    if (!is.character(database))
      stop("`database` must be a character vector of one or more databases ",
           "to test, selected from the following options: ",
           paste(dQuote(choices), collapse = ", "),
           ".")

    database <- match.arg(arg = toupper(database),
                          choices = choices,
                          several.ok = TRUE)

    index <- MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES[database]

  } else {
    # Read molecular signatures from one or more GMT files
    if (!is.character(path_to_gmt))
      stop("`database` must be a character vector of one or more databases ",
           "to test, selected from the following options: ",
           paste(dQuote(choices), collapse = ", "),
           ".")

    index <- lapply(path_to_gmt, TMSig::readGMT)
  }

  names(index) <- NULL # index is a nested list
  index <- unlist(index, recursive = FALSE)

  return(index)
}


#' @title Prepare a list of molecular signatures
#'
#' @param background_list a named list of character vectors specifying the
#'   background used to filter each set in \code{index}.
#' @param index named list of molecular signatures.
#' @inheritParams run_cameraPR
#'
#' @returns An object of class \code{data.frame}.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom data.table data.table :=
#' @importFrom dplyr %>% mutate arrange desc across contains filter n select
#'   pull bind_rows everything where right_join
#' @importFrom tidyr pivot_wider unnest
#'
#' @noRd

.prepare_sets <- function(background_list,
                          index,
                          overlap_cutoff = 0.7,
                          min_size = 5L)
{
  overlap_cutoff <- max(0, min(1, overlap_cutoff))

  # Empty lists to store results
  index_list <- similar_sets_list <- list()

  for (name_i in names(background_list)) {
    ome_i <- sub(".*\\.", "", name_i)

    # Background of genes, metabolites, or phosphosites from the DA results
    background_i <- background_list[[name_i]]

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
      set_dt <- subset(set_dt, subset = elements2 %in%  background_i)

      index_i <- split(x = set_dt[["elements"]], f = set_dt[["sets"]])
      set_sizes <- lengths(index_i)
      keep_sizes <- (set_sizes >= min_size) & (set_sizes < length(background_i))
      index_i <- index_i[keep_sizes]
    } else {
      # Restrict sets to background, filter by size
      index_i <- TMSig::filterSets(
        x = index,
        background = background_i,
        min_size = min_size,
        max_size = length(background_i) - 1L
      )
    }

    # Overlap filter is not used for phosphoproteomics or
    # metabolomics/lipidomics datasets
    if (!ome_i %in% c("prot-ph", "metab")) {
      # Require a minimum proportion of genes in each set to be in the
      # background
      overlap_prop <- lengths(index_i) / lengths(index)[names(index_i)]

      keep <- which(overlap_prop >= ifelse(ome_i == "prot-ol",
                                           0.1,
                                           overlap_cutoff))

      if (length(keep) == 0L) {
        warning("No sets pass `min_size` and `overlap_cutoff` filters for ",
                ome_i,
                ". Excluding this dataset from the results.")

        next
      }

      index_i <- index_i[keep]
    }

    index_list[[name_i]] <- index_i
  }

  return(index_list)
}


#' @title Prepare DA results and molecular signatures
#'
#' @description This function combines .prepare_DA_results and .prepare_sets. It
#'   is used to prepare differential analysis results and molecular signatures
#'   for both CAMERA-PR and ssGSEA.
#'
#' @inheritParams run_cameraPR
#'
#' @returns a named list containing several objects needed for subsequent
#'   analyses.
#'
#' @author Tyler Sagendorf
#'
#' @noRd

.prepare_DA_results_and_sets <- function(DA_list = NULL,
                                         selected_omes = c("transcript-rna-seq",
                                                           "prot-pr",
                                                           "prot-ph",
                                                           "prot-ol",
                                                           "metab"),
                                         selected_tissues = "all",
                                         database = setdiff(
                                           names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
                                           "PTMSIGDB"
                                         ),
                                         path_to_gmt = NULL,
                                         min_size = 5L,
                                         overlap_cutoff = 0.7)
{
  on.exit(gc())

  index <- .create_index(database = database,
                         path_to_gmt = path_to_gmt)

  ## Prepare DA results and molecular signatures ----
  DA_list <- .prepare_DA_results(DA_list = DA_list,
                                 selected_omes = selected_omes,
                                 selected_tissues = selected_tissues,
                                 convert_features = TRUE)

  index_list <- .prepare_sets(background_list = lapply(DA_list, rownames),
                              index = index,
                              min_size = min_size,
                              overlap_cutoff = overlap_cutoff)

  # Collect results in a list
  out <- list("DA_list" = DA_list,
              "index" = index,
              "index_list" = index_list)

  return(out)
}
