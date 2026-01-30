# So the goal of this list is to track the exact type of tissue used in each
# assay. For example, using plasma vs whole blood makes a big difference in the
# interpretation of results or comparison of results across different omes.
# So we make sure that each file has the exact tissue (e.g plasma) and tissue
# code number (e.g. t02) as well as the overarching category (blood, muscle, adipose)

gsutil = "gsutil"
all_gsutil_files = system(paste(gsutil, "ls -R gs://motrpac-data-hub/human-precovid/results"), intern = TRUE)
all_gsutil_files = all_gsutil_files[!base::grepl("t10.*atac-seq|atac-seq.*t10", all_gsutil_files)] #remove t10 atac because that's just the reference standards

tissue_combinations <- list(
  blood = c('t02|t03|t04|t05'),
  muscle = c('t06|t10'),
  adipose = c('t07|t11'))

data_list = lapply(all_gsutil_files, function(path) {
  # Extract tissue code and ome
  tissue_code = sub(".*/(t\\d{2}-[a-zA-Z0-9-]+)/.*", "\\1", path)
  if (tissue_code == path) return(NULL) #basically if the t00 whatever isn't found

  ome_pattern = paste(ome_available_list(), collapse = "|")
  ome = ifelse(base::grepl(ome_pattern, path), base::regmatches(path, base::regexpr(ome_pattern, path)), NA)

  tissue = sapply(names(tissue_combinations), function(name) {
    if (base::grepl(paste(tissue_combinations[[name]], collapse = "|"), tissue_code)) return(name)
    return(NULL)
  })
  tissue = tissue[!sapply(tissue, is.null)][1]; names(tissue) = "tissue"
  if(is.null(tissue_code) | is.null(ome) | is.null(tissue)) return(NULL)
  data.frame("tissue" = tissue, "tissue_code" = tissue_code, "ome" = ome)
})

data_list = data_list[!sapply(data_list, function(x) is.null(x) || any(is.na(x)))]
OME_TISSUE_CODE = do.call(rbind, data_list) %>% dplyr::distinct(., .keep_all = TRUE)
# here we have to manually add the methylseq stuff because its not processed at the BIC

OME_TISSUE_CODE = rbind(OME_TISSUE_CODE, c("blood", "t02-plasma", "metab-meta-reg"))
OME_TISSUE_CODE = rbind(OME_TISSUE_CODE, c("muscle", "t10-muscle", "metab-meta-reg"))
OME_TISSUE_CODE = rbind(OME_TISSUE_CODE, c("adipose", "t11-adipose", "metab-meta-reg"))

OME_TISSUE_CODE = rbind(OME_TISSUE_CODE, c("blood", "t03-edta", "epigen-methylcap-seq"))
OME_TISSUE_CODE = rbind(OME_TISSUE_CODE, c("muscle", "t06-muscle", "epigen-methylcap-seq"))
OME_TISSUE_CODE = rbind(OME_TISSUE_CODE, c("adipose", "t11-adipose", "epigen-methylcap-seq"))

usethis::use_data(OME_TISSUE_CODE, overwrite = TRUE)
