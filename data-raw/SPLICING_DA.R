# SPLICING_DA is an R data object containing results from a differential alternative splicing analysis,
# stored and distributed to facilitate downstream integration and reproducibility.

# The analytical methods used to generate this object—including read alignment, isoform quantification,
# splicing event detection, and statistical testing—are described in detail in
#{"Exercise modulation of the alternative splicing landscape in human tissues"}
# and are not rederived here.

# To accommodate file size and distribution constraints,this object includes only splicing isoforms
# that met the significance threshold FDR < 0.05

# For any questions please reach out to: Zidong Zhang @ zidong.zhang@mssm.edu
SPLICING_DA = readRDS("differential_splicing.rds")

usethis::use_data(SPLICING_DA,
                  overwrite = TRUE)
