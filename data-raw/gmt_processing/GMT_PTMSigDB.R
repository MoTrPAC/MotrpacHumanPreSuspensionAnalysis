library(dplyr)
library(TMSig)

ptmsigdb <- file.path(
  "data-raw",
  "gmt_processing",
  "gmt_files",
  "ptm.sig.db.all.flanking.human.v2.0.0.gmt.gz"
) %>%
  readGMT() %>%
  lapply(function(xi) {
    # Subset to phosphosites only
    xi <- xi[grepl("-p;[ud]$", xi)]

    # Remove phospho indicator
    xi <- sub("-p", "", xi)

    # Replace "_" with "-" and convert middle AA to lower case
    xi <- gsub("_", "-", xi)
    xi <- sub("(^.{7})(.{1})(.*$)", "\\1\\L\\2\\E\\3",
              xi, perl = TRUE)

    return(xi)
  })

names(ptmsigdb) <- paste0("PTMSIGDB_", names(ptmsigdb))

# Write named list to a GMT file
path <- file.path("data-raw", "gmt_processing", "gmt_files", "ptmsigdb.v2.0.flanking.gmt")

MotrpacHumanPreSuspensionAnalysis:::.writeGMT(x = ptmsigdb,
                                      path = file)

# Compress file
R.utils::gzip(file)
