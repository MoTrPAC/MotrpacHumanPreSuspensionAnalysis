library(dplyr)
tissue_mapping = MotrpacHumanPreSuspensionAnalysis::OME_TISSUE_CODE %>%
  select(-ome) %>%
  dplyr::distinct(tissue_code, .keep_all = T)
config = jsonlite::fromJSON("~/config.json")
#basically wherever you want to downloads stuff to.
repo_local_dir = config$precovid_repo_path

metab_file_path = "gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3/resources/motrpac_human-precovid_metabolite-cv_v1.2.txt"
METABOLOMICS_CVS = MotrpacBicQC::dl_read_gcp(path = metab_file_path,
                                      tmpdir = file.path(repo_local_dir, "data", "tmp"))

usethis::use_data(METABOLOMICS_CVS, overwrite = TRUE)
