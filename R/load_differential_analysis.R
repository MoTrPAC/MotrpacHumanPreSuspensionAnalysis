#' @title Load Differential Analysis Results
#'
#' @param selected_omes character; one of \code{\link{ome_available_list}}.
#' @param selected_tissues character; one of
#'   \code{\link{tissue_available_list}}.
#' @param single_matrix logical; if \code{TRUE}, returns a single
#'   \code{data.frame} containing all results. Otherwise, returns a list of
#'   \code{data.frame} objects (default).
#' @param epigen logical; a toggle of TRUE/FALSE if epigenetics data is desired. Loading epigenetic data files is through AWS and is very slow due to file sizes.
#' @param combine_with_featgene logical; whether to include columns from
#'   \code{HUMAN_FEATURE_TO_GENE} in the output.
#' @param verbose logical; whether or not to display messages for some warnings.
#'
#' @returns A nested list of \code{data.table} objects. The top level names are
#'   the tissues, while the second level names are the omes. Each table may
#'   possess the following columns:
#'
#'   \describe{
#'     \item{tissue}{factor; the tissue.}
#'     \item{assay}{factor; the ome.}
#'     \item{platform}{factor; (metabolomics only) metabolomics platform.}
#'     \item{full_model}{factor; full model containing predictors and any
#'     covariates.}
#'     \item{contrast}{factor; full contrast (up to 33).}
#'     \item{contrast_short}{factor; shortened version of the contrasts.}
#'     \item{contrast_type}{factor; one of "exercise_with_controls",
#'     "exercise_no_controls", "Endur_vs_Resist", "baseline", or
#'     "control_only".}
#'     \item{contrast_category}{factor; one of "EE-CON", "RE-CON", "EE-EE",
#'     "RE-RE", "EE-RE", or "CON-CON".}
#'     \item{feature_id}{factor; feature identifier (may be proteins,
#'     phosphosites, transcripts, metabolites/lipids, peaks, or GpGs).}
#'     \item{logFC}{numeric; difference between the group means in the
#'     contrast.}
#'     \item{CI.L}{numeric; lower confidence limit.}
#'     \item{CI.R}{numeric; upper confidence limit.}
#'     \item{degrees_of_freedom}{numeric; degrees of freedom.}
#'     \item{logLik}{numeric; log likelihood of differential expression.}
#'     \item{AveExpr}{numeric; mean of all sample-level values for that
#'     feature.}
#'     \item{methylation_diff}{numeric; methylation difference.}
#'     \item{t}{numeric; moderated t-statistic.}
#'     \item{z.std}{numeric; standard normal equivalent of the t-statistic
#'     (z-scores).}
#'     \item{p_value}{numeric; p-value.}
#'     \item{adj_p_value}{numeric; p-values adjusted within each combination of
#'     tissue, assay, platform, and contrast using the Benjamini-Hochberg method
#'     to control the false discovery rate.}
#'   }
#'
#' @author Tyler Sagendorf Christopher Jin
#'
#' @importFrom data.table setorderv setcolorder rbindlist
#' @importFrom utils data
#'
#' @export load_differential_analysis
#'
#' @examples
#' DA_list <- load_differential_analysis() # default behavior
#'
#' # Structure of a single object
#' str(DA_list[["adipose"]][["prot-pr"]])
#'
#' # Un-nest list
#' DA_list <- unlist(DA_list, recursive = FALSE)
#' names(DA_list)
#'
#' # Include epigen data
#' \dontrun{
#' repo_local_dir <- "path/to/some/directory"
#' DA_list <- load_differential_analysis(
#'   repo_local_dir = repo_local_dir,
#'   epigen = TRUE
#' )
#' }
#'

load_differential_analysis <- function(selected_omes = "all",
                                       selected_tissues = "all",
                                       single_matrix = FALSE,
                                       epigen = FALSE,
                                       combine_with_featgene = FALSE,
                                       verbose = TRUE) {
  selected_tissues <- match.arg(
    arg = selected_tissues,
    choices = c("all", "adipose", "blood", "muscle"),
    several.ok = TRUE
  )

  #-----here I basically just make sure that if any metab platform is listed,
  #all metab is loaded, to support differences in platform specific loading
  if(any(grepl("metab", selected_omes))){
    selected_omes = selected_omes[-grep("metab", selected_omes)]
    selected_omes = c(selected_omes, "metab")
  }

  selected_omes <- match.arg(
    arg = selected_omes,
    choices = c(
      "all", "transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol", "metab",
      "epigen-atac-seq", "epigen-methylcap-seq"
    ),
    several.ok = TRUE
  )

  if ("all" %in% selected_tissues) {
    selected_tissues <- c("adipose", "blood", "muscle")
  }

  # Handle epigen omes requested (explicitly or via "all") while epigen is off
  epigen_omes <- c("epigen-atac-seq", "epigen-methylcap-seq")
  requested_epigen <- intersect(selected_omes, epigen_omes)
  if (!epigen & length(requested_epigen) > 0 &
      all(selected_omes %in% epigen_omes)) {
    stop(
      "You've requested only epigenetic omes (",
      paste(requested_epigen, collapse = ", "),
      ") but `epigen = FALSE`. Set `epigen = TRUE` to load epigenetic data."
    )
  }
  if (verbose & !epigen &
      (length(requested_epigen) > 0 | "all" %in% selected_omes)) {
    message(
      "You've requested one or more epigenetic omes (via explicit selection ",
      "or \"all\") but `epigen = FALSE`, so epigenetic data will be skipped. ",
      "Set `epigen = TRUE` to load epigenetic data."
    )
  }

  if ("all" %in% selected_omes) {
    selected_omes <- c(
      "transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol", "metab",
      "epigen-atac-seq", "epigen-methylcap-seq"
    )
  }

  if (epigen) {
    selected_omes_epigen <- selected_omes[selected_omes %in%
                                            c("epigen-atac-seq",
                                              "epigen-methylcap-seq")]
    if(verbose){
      message("You've elected to load in the epigenetic data too. These file sizes are significantly larger and will require loading in data from AWS. This loading can be quite slow.")
    }
  }

  if(verbose & "metab" %in% selected_omes){
    message("Please remember that the lowest CV Metabolite is chosen and the
            relevant refmet name is used. If you're not able to find your desired
            metabolite, look through the METABOLOMICS_CV object for the relevant
            refmet/feature name.")
  }
  # Split epigen platforms from non-epigen, load epigen via old functionality
  selected_omes <- selected_omes[!selected_omes %in%
                                   c("epigen-atac-seq",
                                     "epigen-methylcap-seq")]

  DA_files <- data(package = "MotrpacHumanPreSuspensionAnalysis")
  DA_files <- DA_files[["results"]][, "Item"]
  DA_files <- DA_files[grepl("_DA$", DA_files)]

  tissues <- tolower(sub("\\_.*", "", DA_files))

  omes <- sub("^[^_]+_(.*)_DA$", "\\1", DA_files)
  omes <- sub("_", "-", tolower(omes))
  omes[omes == "trnscrpt"] <- "transcript-rna-seq"

  new_names <- structure(
    .Data = paste0(tissues, ".", omes),
    names = DA_files
  )

  keep <- (tissues %in% selected_tissues) & (omes %in% selected_omes)
  new_names <- new_names[keep]

  # Load DA results into a list. Only works because of lazy loading
  out <- vector(mode = "list", length = length(new_names))
  names(out) <- as.character(new_names)

  for (i in seq_along(new_names)) {
    out[[i]] <- eval(parse(text = names(new_names[i])))
  }

  if (epigen) {
    epi_list <- load_DA_from_AWS(selected_tissues = selected_tissues,
                                 selected_omes = selected_omes_epigen)

    epi_list <- unlist(epi_list, recursive = FALSE)
    epi_list <- .process_raw_DA(epi_list)
    out <- c(out, epi_list)
  }

  if (combine_with_featgene) {
    out <- lapply(out, function(xi) {
      cols <- colnames(xi)

      xi <- merge(
        x = xi,
        y = MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE,
        by = c("assay", "feature_id"),
        all.x = TRUE,
        all.y = FALSE
      )

      setcolorder(xi, neworder = cols) # reset column order
      setcolorder(
        x = xi,
        neworder = setdiff(
          x = colnames(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE),
          y = cols
        ),
        after = "feature_id"
      )

      # Remove columns with only missing values
      keep_cols <- vapply(xi, function(col_i) any(!is.na(col_i)), logical(1L))

      xi <- xi[, which(keep_cols), with = FALSE]
    })
  }

  if (single_matrix) {
    out <- rbindlist(l = out, use.names = TRUE, fill = TRUE)

    setorderv(x = out, cols = "contrast", order = 1L)
  } else {
    # Nest by tissue
    tissues <- sub("\\..*$", "", names(out))
    names(out) <- sub(".*\\.", "", names(out))
    out <- split(do.call(list, out), tissues)
  }

  return(out)
}


## Internal functions ----------------------------------------------------------

#' @title Process Raw DA Results
#'
#' @description Converts each \code{data.frame} in a list of DA results to a
#'   \code{data.table}, converts character columns to factors, reorders columns,
#'   and sets the key.
#'
#' @param DA_list a named list of DA results (individual \code{data.frame}
#'   objects).
#'
#' @returns A modified version of \code{DA_list} where each list element is a
#'   keyed \code{data.table}. The object will use significantly less memory.
#'
#' @importFrom dplyr %>% left_join select arrange mutate across any_of relocate everything
#' @importFrom data.table as.data.table := setcolorder setorderv setkeyv copy setDT
#'
#' @author Tyler Sagendorf Christopher Jin
#'
#'
#' @noRd

.process_raw_DA <- function(DA_list) {
  for (name_i in names(DA_list)) {
    dt <- copy(DA_list[[name_i]])
    setDT(x = dt)

    # Add contrast information
    dt <- merge(
      x = dt, y = MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER, by = "contrast",
      all.x = TRUE, all.y = FALSE
    )

    contrast_levels <- levels(x = MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER[["contrast"]])

    dt[, `:=`(
      tissue = sub("\\..*$", "", name_i),
      assay = sub("^.*\\.", "", name_i)
    )]

    dt[, `:=`(
      tissue = as.character(tissue),
      assay = as.character(assay),
      contrast = factor(x = contrast, levels = contrast_levels),
      feature_id = as.character(feature_id),
      full_model = as.factor(full_model)
    )]

    dt[, contrast_order := NULL]

    if ("platform" %in% colnames(dt)) {
      dt[, platform := as.factor(platform)]
    }

    dt[, contrast := droplevels(contrast)]

    # Reorder columns
    new_order <- intersect(
      x = c(
        "tissue", "assay", "platform", "full_model",
        colnames(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER)
      ),
      y = colnames(dt)
    )
    setcolorder(x = dt, neworder = new_order)

    if ("z.std" %in% colnames(dt)) {
      setcolorder(x = dt, neworder = "z.std", before = "p_value")
    }

    # Reorder rows
    keys <- intersect(
      x = c(
        "full_model", "contrast", "platform",
        "p_value", "feature_id"
      ),
      y = colnames(dt)
    )
    setkeyv(x = dt, cols = keys)

    DA_list[[name_i]] <- dt
  }

  return(DA_list)
}



#' @title Download Differential Analysis Results from Google Cloud Bucket
#' @description currently not being used. Was previously used for internal consortium members. Not being fully removed because the saving of data objects is still implemented using this function
#' @param repo_local_dir character; path to the local directory. If this
#'   directory does not contain a data/tmp/ subdirectory, one will be created
#'   and files will be downloaded from the appropriate GCP Bucket; otherwise,
#'   files will be read from the directory (or downloaded, if any are missing).
#' @param selected_omes character; one of \code{\link{ome_available_list}}.
#' @param selected_tissues character; one of
#'   \code{\link{tissue_available_list}}.
#' @param single_matrix logical; if \code{TRUE}, returns a single
#'   \code{data.frame} containing all results. Otherwise, returns a list of
#'   \code{data.frame} objects (default).
#' @param epigen logical; whether to download the epigen results. It will take
#'   30 minutes or more.
#' @param gsutil character; the gsutil command.
#' @param load_acute_only logical;
#' @param remove_redundant_metab logical
#' @param include_metab_meta_analysis logical
#' @param combine_with_featgene logical; whether to include columns from
#'   \code{HUMAN_FEATURE_TO_GENE}.
#'
#' @returns A nested list of \code{data.frame} objects or a single
#'   \code{data.frame} containing the differential analysis results.
#'
#' @author Christopher Jin
#'
#' @importFrom data.table rbindlist
#' @importFrom dplyr %>% mutate filter pull arrange left_join
#' @importFrom MotrpacBicQC dl_read_gcp
#' @importFrom stats p.adjust
#'
#' @noRd

.load_differential_analysis <- function(repo_local_dir = NULL,
                                        selected_omes = "all",
                                        selected_tissues = "all",
                                        epigen = FALSE,
                                        gsutil = "gsutil",
                                        load_acute_only = TRUE,
                                        remove_redundant_metab = TRUE,
                                        include_metab_meta_analysis = FALSE)
{

  if (is.null(repo_local_dir)) {
    warning(
      "`repo_local_dir` is not specified, so the current ",
      "working directory will be used.",
      immediate. = TRUE
    )

    repo_local_dir <- getwd()
  }

  tmpdir <- file.path(repo_local_dir, "data", "tmp")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

  da_gsutil_path <- "gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3"

  # Bottleneck 1
  gsutil_files <- system(
    command = paste0(gsutil, " ls -R ", da_gsutil_path),
    intern = TRUE
  )

  gsutil_files <- gsutil_files[grep("*\\.txt$", gsutil_files)]
  gsutil_files <- gsutil_files[grep("_da_", gsutil_files)]

  if (all(selected_omes == "all")) {
    selected_omes <- ome_available_list()
  }

  if (!all(selected_omes %in% ome_available_list())) {
    message(
      "Invalid ome selection. ",
      "Try 'ome_available_list()' for a list of valid omes."
    )
  }

  if (!epigen) {
    selected_omes <- setdiff(
      selected_omes, c("epigen-atac-seq", "epigen-methylcap-seq")
    )
  }

  if (!include_metab_meta_analysis) {
    selected_omes <- setdiff(selected_omes, c("metab-meta-reg"))
  }

  if (all(selected_tissues == "all")) {
    selected_tissues <- tissue_available_list(verbose = FALSE)
  }

  if (!all(selected_tissues %in% tissue_available_list(verbose = FALSE))) {
    message(
      "Invalid tissue selection. ",
      "Try `tissue_available_list()` for a list of valid tissues."
    )
  }

  da_results <- list()

  for (file_path in gsutil_files) {
    tissue <- .find_tissue(file_path)

    if (is.null(tissue) ||
        !tissue %in% selected_tissues) {
      next
    }

    ome <- .find_ome(file_path)

    if (is.null(ome) || !ome %in% selected_omes) {
      next
    }

    # Ensure the structure exists
    if (is.null(da_results[[tissue]])) {
      da_results[[tissue]] <- list()
    }

    if (is.null(da_results[[tissue]][[ome]])) {
      da_results[[tissue]][[ome]] <- list()
    }

    # Bottleneck 2
    file_loaded <- MotrpacBicQC::dl_read_gcp(
      path = file_path,
      tmpdir = file.path(repo_local_dir, "data", "tmp"),
      gsutil = gsutil
    )

    file_loaded$tissue <- tissue

    if (remove_redundant_metab &
        grepl("metab", ome) &
        !grepl("metab-meta-reg", ome)) {
      file_loaded = .prioritize_metab_by_cv_da(file_loaded, tissue, ome)
    }
    da_results[[tissue]][[ome]] <- file_loaded
  }
  return(da_results)
}




#' @title Prioritize redundancies using CVs, Uses refmet names
#'
#' @param file_loaded a data frame of metabolomics differential abundance with
#' the column name structure specified above.
#'
#' @description
#' There are also some situations where multiple feature_ids within one platform
#' map to one given refmet name (e.g. metabolite_isomer1::refmet name & metabolite_isomer2::refmet name).
#' For situations like this, I will be also just choosing lowest CV amongst the
#' options and replacing the name with the refmet name.
#'
#' Every time a load_qc or load_da is called, features should be using refmet
#' names and all spelling/naming inconsistencies will be using refmet stuff.
#'
#' @returns a data frame with no redudant metabolites at a refmet level.
#'
#' @noRd
.prioritize_metab_by_cv_da = function(file_loaded,
                                      tissue,
                                      ome){
  lowest_cv_check = MotrpacHumanPreSuspensionAnalysis::METABOLOMICS_CVS %>%
    dplyr::filter(tissue == !!tissue,
                  assay == ome) %>%
    dplyr::select(feature_id, refmet_name, lowest_CV)
  file_loaded = file_loaded %>%
    dplyr::left_join(., lowest_cv_check, by = "feature_id")  %>%
    dplyr::filter(lowest_CV == "yes") %>%
    dplyr::group_by(contrast) %>%
    dplyr::mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(feature_id = refmet_name) %>%
    dplyr::select(-c(refmet_name, lowest_CV)) %>%
    dplyr::filter(!is.na(feature_id))
  return(file_loaded)
}

