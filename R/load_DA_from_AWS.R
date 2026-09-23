# Public CloudFront release of the c2.0 epigenomics differential analysis.
.AWS_EPIGEN_DA_URL = "https://d1yw74buhe0ts0.cloudfront.net/data/analysis/human_presuspension_sed_adu/c2.0/epigenomics/da/"

# Method tag and file version of each epigenomics DA table in the c2.0 release.
# Versions are per file; the CDN has no directory listing to resolve them from.
.AWS_EPIGEN_DA_FILES = list(
  "epigen-atac-seq" = c(method = "dream-acute", version = "2.1"),
  "epigen-methylcap-seq" = c(method = "malax-glmm-acute", version = "1.2")
)

#' @title Load Epigenomic Differential Analysis Results from AWS
#'
#' @description
#' Loads the publicly released epigenomic differential analysis (DA) results
#' for the selected tissues and omes from the MoTrPAC CloudFront distribution
#' over HTTPS. No credentials are required. Only the epigenomics tier is read
#' this way; every other ome ships inside the package.
#'
#' Files are downloaded on every call and are not cached. The muscle ATAC-seq
#' table is about 2.5 GB.
#'
#' @param selected_tissues character vector; tissues to retrieve.
#' @param selected_omes character vector; epigenomic omes to load (e.g.
#'   \code{"epigen-methylcap-seq"}, \code{"epigen-atac-seq"}).
#'
#' @returns
#' A nested list of \code{data.frame} objects: tissues at the top level, omes at
#' the second. Tissue-ome combinations that were not measured are absent.
#'
#' @author Christopher Jin
#'
#' @seealso \code{\link{load_differential_analysis}}
#'
#' @noRd
.load_DA_from_AWS = function(selected_tissues,
                             selected_omes) {
  epigen_list = list()
  for (tissue in selected_tissues) {
    for (ome in selected_omes) {
      loaded = .load_single_ome_tissue_AWS(selected_ome = ome,
                                           selected_tissue = tissue)
      if (is.null(loaded)) next
      epigen_list[[tissue]][[ome]] = loaded
    }
  }
  return(epigen_list)
}

#' @title URL of a Single Epigenomic DA Table on AWS
#'
#' @param selected_ome character; epigenomic ome.
#' @param selected_tissue character; tissue identifier.
#'
#' @returns A length-1 character URL, or \code{NULL} if the tissue-ome
#'   combination was not measured.
#'
#' @noRd
.aws_epigen_da_url = function(selected_ome,
                              selected_tissue) {
  tissue_code = MotrpacHumanPreSuspensionAnalysis::OME_TISSUE_CODE %>%
    dplyr::filter(ome == selected_ome,
                  tissue == selected_tissue) %>%
    dplyr::pull(tissue_code)
  if (length(tissue_code) == 0) return(NULL)

  file_spec = .AWS_EPIGEN_DA_FILES[[selected_ome]]
  file_name = paste("human-precovid-sed-adu", tissue_code, selected_ome, "da",
                    file_spec[["method"]], sep = "_")
  file_url = paste0(.AWS_EPIGEN_DA_URL, file_name, "_v", file_spec[["version"]],
                    ".txt")
  return(file_url)
}

#' @title Load a Single Epigenomic DA Table from AWS
#'
#' @param selected_ome character; epigenomic ome to load.
#' @param selected_tissue character; tissue identifier.
#'
#' @returns A \code{data.frame}, or \code{NULL} if the tissue-ome combination
#'   was not measured.
#'
#' @importFrom utils read.csv
#'
#' @noRd
.load_single_ome_tissue_AWS = function(selected_ome,
                                       selected_tissue) {
  file_url = .aws_epigen_da_url(selected_ome = selected_ome,
                                selected_tissue = selected_tissue)
  if (is.null(file_url)) return(NULL)

  old_timeout = options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old_timeout), add = TRUE)

  loaded_file = read.csv(file_url,
                         sep = "\t",
                         check.names = FALSE)
  return(loaded_file)
}
