library(dplyr)
# File downloaded 2025-01-23 from https://maayanlab.cloud/Enrichr/#libraries and
# compressed with R.utils::gzip
path <- file.path("data-raw", "gmt_processing", "CellMarker_2024.txt.gz")

cellmarker <- TMSig::readGMT(path, check = FALSE)
cellmarker <- cellmarker[grepl("human", names(cellmarker),
                               ignore.case = TRUE)]

names(cellmarker) <- paste0("CELLMARKER_", names(cellmarker))

# Write named list to a GMT file
gmt_file <- file.path("data-raw", "gmt_processing", "gmt_files",
                      "cellmarker.v2024.symbols.gmt")

MotrpacHumanPreSuspensionAnalysis:::.writeGMT(x = cellmarker,
                                            path = gmt_file)

# Compress file
R.utils::gzip(gmt_file)
