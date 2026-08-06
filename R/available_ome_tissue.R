#' @title List All Omes
#'
#' @description Return a character vector containing all available omes. These
#'   are the "full names" for each of the assays.
#'
#' @details This vector is the vocabulary every accessor gates on:
#'   \code{load_qc()}, \code{load_summary_stats()} and
#'   \code{load_differential_analysis()} all resolve \code{"all"} through it and
#'   reject anything absent from it, so an ome missing here is an ome no
#'   accessor can return. It must therefore track what the pipeline actually
#'   builds, which is recorded in \code{OME_TISSUE_CODE}.
#'
#'   Two clinical omes are included as of v2.0. Clinical chemistry was a single
#'   \code{"clinical-chemistry"} assay in v1.3 and is now split into
#'   \code{"metab-t-clinical"} and \code{"prot-clinical"}, each with its own QC,
#'   differential-analysis and summary-statistic objects.
#'
#'   \code{"metab-meta-reg"} was removed. It was dropped from the pipeline
#'   entirely, so no object was ever produced under that name and every accessor
#'   offered it as a choice that returned nothing.
#'
#'   The \code{lab-*} tiers in \code{OME_TISSUE_CODE} are deliberately absent.
#'   They are the raw clinical laboratory inputs, distributed as the
#'   \code{cln_chemistry_*} tables rather than as omes.
#'
#' @returns A character vector containing all available omes.
#'
#' @seealso [metab_only_list()], [tissue_available_list()]
#'
#' @export ome_available_list
#'
#' @examples
#' ome_available_list()
ome_available_list <- function() {
  out <- c(
    "prot-ol", "prot-ph", "prot-pr", "prot-clinical", "transcript-rna-seq",
    "epigen-methylcap-seq", "epigen-atac-seq",
    "metab-u-hilicpos", "metab-u-ionpneg", "metab-u-lrpneg", "metab-u-lrppos",
    "metab-u-rpneg", "metab-u-rppos", "metab-t-amines", "metab-t-conv",
    "metab-t-imm-crt", "metab-t-oxylipneg", "metab-t-tca", "metab-t-nuc",
    "metab-t-acoa", "metab-t-ka", "metab-t-clinical"
  )

  return(out)
}


#' @title List the Metabolomics Platforms
#'
#' @description A character vector containing the targeted and untargeted
#'   metabolomics platforms.
#'
#' @details These are the research metabolomics platforms.
#'
#'   \code{"metab-t-clinical"} is deliberately NOT here, though it is a
#'   targeted metabolomics platform and a valid choice in
#'   \code{\link{ome_available_list}()}. It is clinical chemistry, and the
#'   loaders gate it behind \code{load_clinical} rather than folding it into
#'   the research platforms — see \code{\link{clinical_ome_list}()}.
#'
#' @returns A character vector containing the targeted and untargeted
#'   metabolomics platforms.
#'
#' @seealso [ome_available_list()], [clinical_ome_list()]
#'
#' @export metab_only_list
#'
#' @examples
#' metab_only_list()
metab_only_list <- function() {
  out <- c(
    "metab-u-hilicpos", "metab-u-ionpneg", "metab-u-lrpneg",
    "metab-u-lrppos", "metab-u-rpneg", "metab-u-rppos",
    "metab-t-amines", "metab-t-conv", "metab-t-imm-crt",
    "metab-t-oxylipneg", "metab-t-tca", "metab-t-nuc", "metab-t-acoa",
    "metab-t-ka"
  )

  return(out)
}

#' @title List the Clinical Omes
#'
#' @description The omes carrying clinical chemistry rather than a research
#'   assay. Every loader takes a \code{load_clinical} argument that gates these,
#'   and it is \code{FALSE} by default.
#'
#' @details v1.3 carried clinical chemistry as a single \code{"clinical-chemistry"}
#'   assay; v2.0 splits it into a metabolomics and a proteomics assay, each with
#'   its own QC, differential-analysis and summary-statistic objects.
#'
#'   They are gated rather than simply included because they are a different
#'   kind of measurement from the research omes, and every caller written before
#'   v2.0 that asks for \code{"all"} is summarising the molecular landscape.
#'   Adding them by default changes those results silently: the clinical
#'   metabolomics differential-analysis rows overlap the combined
#'   \code{*_METAB_DA} table on five analytes (Cortisol, Glycerol, KET, NEFA and
#'   Glucose), and a caller that maps an assay to a display name puts clinical
#'   chemistry into a row that was never meant to hold it.
#'
#'   Pass \code{load_clinical = TRUE} to get them.
#'
#' @returns A character vector of the clinical omes.
#'
#' @seealso [ome_available_list()], [metab_only_list()]
#'
#' @export clinical_ome_list
#'
#' @examples
#' clinical_ome_list()
clinical_ome_list <- function() {
  out <- c("prot-clinical", "metab-t-clinical")

  return(out)
}


#' @title List the Available Tissues
#'
#' @description Vector of tissues available in the analysis. These are the "full
#'   names" for each of the tissues.
#'
#' @param verbose logical; whether to output messages.
#'
#' @returns A character vector of the available tissues: "adipose", "blood", and
#'   "muscle".
#'
#' @export tissue_available_list
#'
#' @examples
#' tissue_available_list()
tissue_available_list <- function(verbose = TRUE) {
  if (verbose)
    message(
      "The available tissues are placed into overarching categories. ",
      "For example, an assay using PBMCs would be categorized as blood."
    )

  out <- c("adipose", "blood", "muscle")

  return(out)
}


#' @title Determine the Ome from a File Path
#'
#' @description Find the corresponding ome within a file path. Helpful for the
#'   gsutil structure we have. Used to split things like the qc-objects into
#'   corresponding ome folders.
#'
#' @param file_path character; path to a file.
#'
#' @returns Character string; the ome contained in the file path, or \code{NULL}
#'   if no match.
#'
#' @examples
#' # Returns "transcript-rna-seq"
#' .find_ome(file.path(
#'   "gs://motrpac-data-hub/human-precovid/results",
#'   "t02-transcript-rna-seq...something"
#' ))
#'
#' @noRd

.find_ome <- function(file_path) {
  omes <- ome_available_list()

  for (ome in omes) {
    if (grepl(ome, file_path)) {
      return(ome)
    }
  }

  # message("No ome found in the file path.")

  return(NULL)
}



#' @title Determine the Tissue from a File Path
#'
#' @description Find the corresponding tissue within a file path. Helpful for
#'   the gsutil structure we have. Used to split things like the qc-objects into
#'   corresponding tissue folders.
#'
#' @param file_path character; path to a file.
#'
#' @returns Character string; the tissue contained in the file path. One of
#'   "adipose", "blood", "muscle", or \code{NULL} if no match.
#'
#' @examples
#' # Returns "blood"
#' .find_tissue(file.path(
#'   "gs://motrpac-data-hub/human-precovid/results",
#'   "t02-transcript...something"
#' ))
#'
#' @noRd

.find_tissue <- function(file_path) {
  tissue_combinations <- list(
    "t02" = "blood",
    "t03" = "blood",
    "t04" = "blood",
    "t05" = "blood",
    "t06" = "muscle",
    "t10" = "muscle",
    "t07" = "adipose",
    "t11" = "adipose"
  )

  for (code in names(tissue_combinations)) {
    if (grepl(code, file_path)) {
      return(tissue_combinations[[code]])
    }
  }

  # message("No tissue code found in the file path.")

  return(NULL)
}

