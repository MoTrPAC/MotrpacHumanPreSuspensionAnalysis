library(dplyr)
library(readxl)

# Human MitoCarta3.0 database
mitocarta <- file.path("data-raw", "gmt_processing", "Human.MitoCarta3.0.xls") %>%
  readxl::read_xls(sheet = "C MitoPathways") %>%
  select(MitoPathway, Genes) %>%
  filter(!is.na(MitoPathway)) %>%
  mutate(Genes = strsplit(Genes, split = ", ")) %>%
  {structure(.$Genes, names = .$MitoPathway)}

names(mitocarta) <- paste0("MITOCARTA_", names(mitocarta))

# Write named list to a GMT file
path <- file.path("data-raw", "gmt_processing", "gmt_files", "mitocarta3.0.symbols.gmt")

MotrpacHumanPreSuspensionAnalysis:::.writeGMT(x = mitocarta,
                                      path = path)

# Compress file
R.utils::gzip(file)
