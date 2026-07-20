library(MotrpacBicQC)
assay_codes = MotrpacBicQC::assay_codes


usethis::use_data(assay_codes,
                  internal = FALSE,
                  overwrite = TRUE,
                  compress = TRUE,
                  version = 3)
