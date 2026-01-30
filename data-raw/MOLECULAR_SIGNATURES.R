library(TMSig) # readGMT

gmt_files <- c(
  # Gene sets from the C2 collection of MSigDB
  "C2" = "c2.all.v2023.2.Hs.symbols.gmt.gz",
  # Gene sets from the C5 Gene Ontology subcollection of MSigDB
  "GO" = "c5.go.v2023.2.Hs.symbols.gmt.gz",
  # MitoCarta3.0 database
  "MITOCARTA" = "mitocarta3.0.symbols.gmt.gz",
  # Kinase sets
  "PSP" = "phosphositeplus.v6.7.1.1.flanking.gmt.gz",
  # RefMet chemical subclasses (metabolomics/lipidomics)
  "REFMET" = "metabolomics.workbench.refmet.2024.08.07.metabolites.gmt.gz",
  # PTM Signatures Database (PTMSigDB)
  "PTMSIGDB" = "ptmsigdb.v2.0.flanking.gmt.gz",
  # CellMarker gene sets
  "CELLMARKER" = "cellmarker.v2024.symbols.gmt.gz"
)

gmt_files <- structure(file.path("data-raw", "gmt_files", gmt_files),
                       names = names(gmt_files))

# Nested list of molecular signatures
MOLECULAR_SIGNATURES <- lapply(gmt_files, readGMT)

# Separate Gene Ontology databases
GO <- MOLECULAR_SIGNATURES[["GO"]]
ont <- sub("^(GO[^_]+)_.*", "\\1", names(GO))
GO <- split(do.call(list, GO), f = ont)

C2 <- MOLECULAR_SIGNATURES[["C2"]]
databases <- c("BIOCARTA", "KEGG_MEDICUS", "PID", "REACTOME", "WP")

keep <- grepl(paste(paste0("^", databases, "_"), collapse = "|"),
              names(C2))

C2 <- C2[keep]
group <- sub("^([^_]+).*", "\\1", names(C2))
group[group == "KEGG"] <- "KEGG_MEDICUS"

C2 <- split(do.call(list, C2), f = group)

MOLECULAR_SIGNATURES[c("C2", "GO")] <- NULL
MOLECULAR_SIGNATURES <- c(C2, GO, MOLECULAR_SIGNATURES)

# Number of gene sets in each database
lengths(MOLECULAR_SIGNATURES)
# BIOCARTA  KEGG_MEDICUS  PID  REACTOME   WP   GOBP  GOCC  GOMF
#      292           619  196      1692  791   7647  1015  1799
# MITOCARTA  PSP  REFMET  PTMSIGDB  CELLMARKER
#       149  412     162       954        1134

# Save
usethis::use_data(MOLECULAR_SIGNATURES, overwrite = TRUE,
                  version = 3, compress = TRUE)
