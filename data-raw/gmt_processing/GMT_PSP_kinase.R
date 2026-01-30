library(dplyr)

# NOTICE: PSP data is not for commercial use.

# Kinase_Substrate_Dataset obtained from PhosphoSitePlus v6.7.1.1
# https://www.phosphosite.org/staticDownloads.action
# Last Modified Fri Nov 17 08:50:20 EST 2023
kinase_sets <- file.path("data-raw",
                         "gmt_processing",
                         "Kinase_Substrate_Dataset.gz") %>%
  gzfile() %>%
  read.delim(skip = 2) %>%
  filter(KIN_ORGANISM == "human", KIN_ORGANISM == SUB_ORGANISM) %>%
  mutate(human_flanking = toupper(`SITE_...7_AA`),
         human_flanking = gsub("_", "-", human_flanking),
         human_flanking = sub("(^.{7})(.{1})(.*$)",
                              "\\1\\L\\2\\E\\3",
                              human_flanking,
                              perl = TRUE)) %>%
  select(human_flanking, kinase = GENE) %>%
  # mutate(human_site = paste(SUB_ACC_ID, SUB_MOD_RSD, sep = "_")) %>%
  # select(human_site, kinase = GENE) %>%
  distinct() %>%
  unstack()

names(kinase_sets) <- paste0("PSP_", names(kinase_sets))

# Write named list to a GMT file
path <- file.path("data-raw", "gmt_processing", "gmt_files",
                  "phosphositeplus.v6.7.1.1.flanking.gmt")

MotrpacHumanPreSuspensionAnalysis:::.writeGMT(x = kinase_sets,
                                      path = path)

# Compress file
R.utils::gzip(file)
