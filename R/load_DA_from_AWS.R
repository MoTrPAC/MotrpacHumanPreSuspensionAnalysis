#' @title Load Epigenomic Differential Analysis Results from AWS
#'
#' @description
#' Loads publicly released epigenomic differential analysis (DA) results for
#' selected tissues and omes directly from an AWS-backed content delivery
#' network (CDN). These datasets are hosted on Amazon CloudFront and are
#' accessed via HTTPS without requiring authentication or AWS credentials.
#' Only needed for epigenetic files because these are too large for the package.
#'
#' This function iterates over the requested tissue–ome combinations and
#' retrieves each corresponding DA table using
#' \code{\link{.load_single_ome_tissue_AWS}}. The returned object mirrors the
#' nested list structure used elsewhere in the package, with tissues at the
#' top level and omes at the second level.
#'
#' The AWS-hosted files correspond to the public release of the human
#' pre-suspension sedentary adult epigenomics analyses.
#'
#' @param selected_omes character vector; one or more epigenomic omes to load
#'   (e.g., \code{"epigen-methylcap-seq"}, \code{"epigen-atac-seq"}).
#' @param selected_tissues character vector; one or more tissues for which
#'   differential analysis results should be retrieved.
#'
#' @returns
#' A nested list of \code{data.frame} objects. The top-level names correspond to
#' tissues, and the second-level names correspond to omes. Each entry contains
#' a differential analysis table loaded from the AWS public release.
#'
#' @details
#' Files are served via Amazon CloudFront and downloaded on demand at runtime.
#' Connection timeouts are increased internally to accommodate large epigenomic
#' result files.
#'
#' @author Christopher Jin
#'
#' @seealso
#' \code{\link{.load_single_ome_tissue_AWS}}
#' \code{\link{load_differential_analysis}}
#'
#' @examples
#' \dontrun{
#' DA_epigen <- load_DA_from_AWS(
#'   selected_omes = c("epigen-methylcap-seq"),
#'   selected_tissues = c("muscle", "adipose")
#' )
#' }
load_DA_from_AWS = function(selected_tissues,
                            selected_omes){
  epigen_list = list()
  for(tissue in selected_tissues){
    for(ome in selected_omes){
      epigen_list[[tissue]][[ome]] = .load_single_ome_tissue_AWS(
        selected_ome = ome,
        selected_tissue = tissue
      )
    }
  }
  return(epigen_list)
}

#' @title Load a Single Epigenomic DA Table from AWS
#'
#' @description
#' Downloads and loads a single epigenomic differential analysis (DA) result
#' table for a specified tissue and ome from the public AWS release. Files are
#' hosted on Amazon CloudFront and accessed via a static HTTPS URL.
#'
#' This function is intended for internal use and is called by
#' \code{\link{load_DA_from_AWS}} to populate tissue–ome combinations.
#'
#' @param selected_ome character; epigenomic ome to load (e.g.,
#'   \code{"epigen-methylcap-seq"} or \code{"epigen-atac-seq"}).
#' @param selected_tissue character; tissue identifier.
#' @param data_category character; analysis category, default is \code{"da"}
#'   for differential analysis.
#' @param version character; dataset version identifier appended to the file
#'   name
#'
#' @returns
#' A \code{data.frame} containing the differential analysis results for the
#' specified tissue and ome. If no matching tissue–ome combination is found,
#' the function returns \code{NULL}.
#'
#' @details
#' The underlying files are part of the public AWS release of the human
#' pre-suspension sedentary adult epigenomics analyses. File naming conventions
#' vary slightly across epigenomic assays and are handled internally.
#'
#' To accommodate large epigenomic result files, the global R connection
#' timeout is temporarily increased within this function.
#'
#' @author Christopher Jin
#'
#' @keywords internal

.load_single_ome_tissue_AWS = function(selected_ome,
                                       selected_tissue,
                                       data_category = "da",
                                       version = "1.2"){
  #the download for the larger epigen files takes quite a while.
  options(timeout = 1200)

  AWS_header = "https://d1yw74buhe0ts0.cloudfront.net/data/analysis/human_presuspension_sed_adu/v1.3/epigenomics/da/"
  all_file_header = "human-precovid-sed-adu" #this is the base structure for all files within the phase.
  tissue_code = MotrpacHumanPreSuspensionAnalysis::OME_TISSUE_CODE %>%
    dplyr::filter(ome == selected_ome,
                  tissue == selected_tissue) %>%
    dplyr::pull(tissue_code)
  if(length(tissue_code) == 0) return(NULL)

  #we had some special naming conventions based on the methods. this is a bit hard coded but since only the epigen data needs to be private,
  #it doesn't end up being too bad.
  if(selected_ome == "epigen-methylcap-seq") data_details = "malax-glmm-acute"
  if(selected_ome == "epigen-atac-seq") data_details = "dream-acute"

  file_name = paste(all_file_header, tissue_code, selected_ome, data_category, data_details, sep = "_")
  file_name = paste0(AWS_header, file_name, "_v", version, ".txt")

  loaded_file = read.csv(file_name,
                         sep = "\t")
  return(loaded_file)
}


