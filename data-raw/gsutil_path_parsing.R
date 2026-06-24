#' @title Internal function to write files according to BIC naming guidelines
#' The goal of this is to avoid hard coding file paths for outputs
#'
#' @param local_path path to actually output the file structure (usually the data folder in the precovid-analyses repo)
#' @param ome The ome (specific--e.g. prot-pr, prot-ph, etc.) Use. \code{ome_available_list} for reference
#' @param data_category The overarching category for the data (only options: DA, qc-norm, metadata)
#' @param data_details the subcategory for the data (options: samples (for metadata_samples), features (for metadata_features))
#' @param actual_data_object the actual object that you plan on writing out (e.g. a data.frame or something)
#' @param tissue Desired tissue in the file path one of\code{\link{tissue_available_list}}
#' @param version the version that's going to be outputted
#' @param return_name_only if you don't plan on actually writing out a file, and
#' only want the file name that satisfies naming conventions.
#'
#' @returns If `return_name_only` is TRUE, this returns the name of the full
#' naming scheme under motrpac precovid standards. If it's F, it just writes
#' the `actual_data_object` to the local_path specified using the naming scheme.
#'
#' @importFrom utils write.table
#'
#' @noRd

write_with_path_name = function(actual_data_object = NULL,
                                local_path = NULL,
                                ome = NULL,
                                tissue = NULL,
                                data_category = NULL,
                                data_details = NULL,
                                version = "1.2",
                                return_name_only = FALSE){
  file_type = ".txt"
  all_file_header = "human-precovid-sed-adu" #this is the base structure for all files within the phase.
  tissue_code = .match_ome_tissue_code(desired_ome = ome, input_tissue = tissue)
  file_name = paste(all_file_header, tissue_code, ome, data_category, data_details, sep = "_")
  file_name = paste0(local_path,"/", file_name, "_v", version, file_type)
  if(return_name_only) {
    return(file_name)
  }else{
    write.table(actual_data_object, file = file_name, row.names = F, sep = '\t', quote = F)
  }
}

#' @title Internal function to find files according to specific naming
#' @description
#' Pulls specific file paths from the results section of the data hub for human
#' precovid according to specific file names.
#'
#' Actual order -> start with this -> check for most recent version -> check for existing files within tmp folder
#' Then prompt for individuals
#'
#' @param tissue Desired tissue in the file path one of\code{\link{tissue_available_list}}
#' @param data_type Specific file types in the file path (e.g. counts)
#' @param version If a specific file version is desired
#' @param desired_ome Desired tissue in the file path; one of\code{\link{ome_available_list}}
#' @param gsutil character; path to the gsutil executable. Defaults to "gsutil",
#'   which assumes the path to the gsutil executable has been added to the PATH
#'   variable.
#'
#' @returns a character vector of file paths from the results section that satisfy
#' the requirements from the function.
#'
#' @examples
#' file_load = .find_path_name(desired_ome, tissue = tissue, data_type = 'count', version = "1.0") #find version 1 of the counts file
#' @noRd
.find_path_name = function(desired_ome,
                           tissue,
                           data_type,
                           version = 'any',
                           gsutil = "gsutil"){
  all_gsutil_files = system(paste0(gsutil, " ls -R gs://motrpac-data-hub/human-precovid/results"), intern = T)
  all_gsutil_files = all_gsutil_files[grep("\\.txt$|\\.csv$|\\.txt.gz", all_gsutil_files)] #first filter for ending files
  all_gsutil_files = all_gsutil_files[grep(desired_ome , all_gsutil_files)]
  if (length(all_gsutil_files) == 0) {message("No files found matching the desired ome criteria.")}

  actual_tissue_code = .match_ome_tissue_code(desired_ome, tissue) #so here sometimes the tissue we generally use in conversation (muscle, adipose) isn't the
  #exact tissue code e.g. (t10-muscle-powder or whatever, this finds the actual tissue code)
  all_gsutil_files = all_gsutil_files[grep(actual_tissue_code , all_gsutil_files)]
  if (length(all_gsutil_files) == 0) {message("No files found matching the desired ome + tissue criteria.")}
  all_gsutil_files = all_gsutil_files[grep(data_type , all_gsutil_files)]
  if (length(all_gsutil_files) == 0) {message("No files found matching the desired ome, tissue, datatype criteria.")}

  if (version == 'any'){return(all_gsutil_files)}
  #continue filtering if you want a specific version number
  all_gsutil_files = all_gsutil_files[grep(version , all_gsutil_files)]
  if (length(all_gsutil_files) == 0) {message("No files found matching the desired ome, tissue, datatype, version criteria.")}
  return(all_gsutil_files)
}


#' An internal function to find the character code (t01/t02...) for a given ome and tissue.
#'
#' @param desired_ome Desired tissue in the file path; one of\code{\link{ome_available_list}}
#' @param tissue Desired tissue in the file path one of\code{\link{tissue_available_list}}
#' @param gsutil character; path to the gsutil executable. Defaults to "gsutil",
#'   which assumes the path to the gsutil executable has been added to the PATH
#'   variable.
#'
#' @returns a character value of the combined tissue code and tissue value
#'
#' @examples
#' .match_ome_tissue_code("transcript-rna-seq", "muscle") #returns t06-muscle
#' @noRd
.match_ome_tissue_code = function(desired_ome, input_tissue){
  specific_output = MotrpacHumanPreSuspensionAnalysis::OME_TISSUE_CODE %>%
    dplyr::filter(ome == desired_ome, tissue == input_tissue)
  if(nrow(specific_output) == 0) stop("No ome matches the desired ome/tissue combo")
  tissue_name_output = specific_output[["tissue_code"]]
  return(tissue_name_output)
}


.validate_staging_folder = function(){
  current_staging_folder = "gs://pre-cawg/staging_20260428"
  system(paste("gsutil cp -R gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3/epigenomics", current_staging_folder))
  system(paste("gsutil cp -R gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3/metabolomics-targeted", current_staging_folder))
  system(paste("gsutil cp -R gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3/metabolomics-untargeted", current_staging_folder))
  system(paste("gsutil cp -R gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3/proteomics", current_staging_folder))
  system(paste("gsutil cp -R gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3/resources", current_staging_folder))
  system(paste("gsutil cp -R gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3/transcriptomics", current_staging_folder))
  system(paste("gsutil cp -R gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3/clinical_chemistry", current_staging_folder))


}




