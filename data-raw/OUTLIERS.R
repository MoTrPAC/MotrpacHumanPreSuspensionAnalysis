gsutil = "gsutil"
config = jsonlite::fromJSON("~/config.json")
repo_local_dir = config$precovid_repo_path

available_files = system(paste0(gsutil, " ls -R gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3"), intern = TRUE)
outlier_files = grep("removed-samples", available_files, value = TRUE)

outlier_list = lapply(outlier_files, function(path) {
  outlier_file_single = MotrpacBicQC::dl_read_gcp(path, sep = '\t', tmpdir = paste0(repo_local_dir, "/data/tmp/"))
  tissue = .find_tissue(path)
  ome = .find_ome(path)

  #----needed because of the different naming conventions in different
  if("sample" %in% names(outlier_file_single)) outlier_file_single = outlier_file_single %>% dplyr::rename(vialLabel = sample)
  if("Sample" %in% names(outlier_file_single)) outlier_file_single = outlier_file_single %>% dplyr::rename(vialLabel = Sample)
  if("PC" %in% names(outlier_file_single)) outlier_file_single = outlier_file_single %>% dplyr::rename(reason = PC)

  outlier_file_single %>% dplyr::select(vialLabel, reason) %>%
    dplyr::mutate(tissue = tissue,
           ome = ome)
})
OUTLIERS = do.call(rbind, outlier_list)
#here we manually annotate a outlier that had adjusted values.
OUTLIERS = rbind(OUTLIERS, c("11263010401", "manual check", "blood", "transcript-rna-seq")) #rna-seq blood
usethis::use_data(OUTLIERS, overwrite = TRUE)


