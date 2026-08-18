#' @title Load Epigenomic Differential Analysis Results from Google Cloud Storage
#'
#' @description
#' Loads epigenomic differential analysis (DA) results for the selected tissues
#' and omes from a Google Cloud Storage prefix. Only the epigenomics tier is
#' read this way; every other ome ships inside the package. The files are far
#' too large to distribute, and the bucket they live in is consortium-gated, so
#' this requires a working \code{gsutil} and read access.
#'
#' The prefix is listed once with \code{gsutil ls -R} and filtered down to
#' differential analysis tables, then each match is resolved to a tissue and ome
#' from its file name. Files are cached under \code{repo_local_dir}, so a second
#' call in the same project reuses the download instead of repeating it.
#'
#' @param selected_tissues character vector; tissues to retrieve.
#' @param selected_omes character vector; epigenomic omes to load (e.g.
#'   \code{"epigen-methylcap-seq"}, \code{"epigen-atac-seq"}).
#' @param repo_local_dir character; local directory used as the download cache.
#'   Files are written under its \code{data/tmp/} subdirectory, which is created
#'   if absent.
#' @param gsutil character; path to the gsutil executable.
#' @param bucket character; the GCS prefix to read.
#'
#' @returns
#' A nested list of \code{data.frame} objects: tissues at the top level, omes at
#' the second. Tissue-ome combinations with no file in the bucket are absent
#' rather than \code{NULL}.
#'
#' @author Christopher Jin
#'
#' @seealso \code{\link{load_differential_analysis}}
#'
#' @keywords internal

.load_DA_from_bucket <- function(selected_tissues,
                                 selected_omes,
                                 repo_local_dir = NULL,
                                 gsutil = "gsutil",
                                 bucket = .STAGING_BUCKET) {
  if (is.null(repo_local_dir)) {
    stop("`repo_local_dir` has to be specified if you're loading epigenetics ",
         "data. It is the local cache the files are downloaded into.")
  }

  # dl_read_gcp is the only thing this package needs MotrpacBicQC for, and it is
  # only reachable behind `epigen = TRUE`. Keeping it in Suggests means the
  # package installs without it.
  if (!requireNamespace("MotrpacBicQC", quietly = TRUE)) {
    stop("Loading epigenetics data requires the MotrpacBicQC package. ",
         "Install it with ",
         "`remotes::install_github(\"MoTrPAC/MotrpacBicQC\")`.")
  }

  tmpdir <- file.path(repo_local_dir, "data", "tmp")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

  message("Epigenomics differential analysis is not shipped in the package. ",
          "It relies on the gsutil implementation: listing '", bucket,
          "' via `", gsutil, " ls -R` and downloading through ",
          "MotrpacBicQC::dl_read_gcp(). This requires gsutil on your PATH and ",
          "consortium read access to the bucket. Files are cached in ", tmpdir,
          ".")

  bucket_files <- system(paste0(gsutil, " ls -R ", bucket), intern = TRUE)
  bucket_files <- bucket_files[grep("\\.txt$", bucket_files)]
  bucket_files <- bucket_files[grep("_da_", bucket_files)]

  da_results <- list()

  for (file_path in bucket_files) {
    tissue <- .find_tissue(file_path)
    if (is.null(tissue) || !tissue %in% selected_tissues) next

    ome <- .find_ome(file_path)
    if (is.null(ome) || !ome %in% selected_omes) next

    if (is.null(da_results[[tissue]])) da_results[[tissue]] <- list()

    da_results[[tissue]][[ome]] <- MotrpacBicQC::dl_read_gcp(
      path = file_path,
      tmpdir = tmpdir,
      gsutil = gsutil
    )
  }

  if (length(da_results) < 1) {
    message("No epigenetic ome/tissue combination was found in ", bucket,
            ". Please double check your input parameters and bucket access.")
  }

  return(da_results)
}
