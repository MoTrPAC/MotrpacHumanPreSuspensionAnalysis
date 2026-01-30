#from: https://humantfs.ccbr.utoronto.ca/download.php -> see R object
UTORONTO_TFs = read.csv("~/Downloads/DatabaseExtract_v_1.01.txt",sep = "\t", check.names = T)[-1] %>%
  rename(gene_symbol = HGNC.symbol) %>%
  dplyr::left_join(., HUMAN_FEATURE_TO_GENE, by = "gene_symbol") %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  filter(platform == "prot-ph") %>%
  filter(Is.TF. == "Yes") %>%
  select(feature_id, gene_symbol)
usethis::use_data(UTORONTO_TFs, overwrite = TRUE)
