## This script will be changed in the future. For now, I am just modifying the
## HUMAN_FEATURE_TO_GENE table from MotrpacHumanPreSuspension.

library(dplyr)
library(data.table)


# # Add flanking sequence to HUMAN_FEATURE_TO_GENE ----
#
# # Load config file
# config <- jsonlite::read_json("~/configs/config-nmclark2.json")
#
# # Load the table
# HUMAN_FEATURE_TO_GENE <- HUMAN_FEATURE_TO_GENE %>%
#   ungroup() %>%
#   rename(assay = platform) %>%
#   mutate(assay = ifelse(assay == "metabolomics", "metab", assay),
#          assay = as.factor(assay))
#
# # Obtain the flanking sequences from QC norm data
# qc <- load_qc(repo_local_dir = config$gitdir,
#               selected_omes = "prot-ph",
#               gsutil = config$gsutil)
#
# adipose.flanking <-
#   qc$adipose$`prot-ph`$feature_metadata[, c("id", "flanking_sequence")]
#
# muscle.flanking <-
#   qc$muscle$`prot-ph`$feature_metadata[, c("id", "flanking_sequence")]
#
# all.flanking <- rbind(adipose.flanking, muscle.flanking) %>%
#   unique() %>%
#   rename(feature_id = id)
#
# # Add flanking sequences
# HUMAN_FEATURE_TO_GENE <- left_join(HUMAN_FEATURE_TO_GENE,
#                                    all.flanking,
#                                    by = "feature_id")

HUMAN_FEATURE_TO_GENE <- MotrpacHumanPreSuspension::HUMAN_FEATURE_TO_GENE

setDT(HUMAN_FEATURE_TO_GENE)

cols <- colnames(HUMAN_FEATURE_TO_GENE)
char_cols <- vapply(cols, class, character(1L)) == "character"
char_cols <- setdiff(names(char_cols),
                     c("relationship_to_gene", "assay"))
HUMAN_FEATURE_TO_GENE[, (char_cols) := lapply(.SD, as.factor),
                      .SDcols = char_cols]

setcolorder(x = HUMAN_FEATURE_TO_GENE,
            neworder = c("assay", "feature_id"))

setkeyv(x = HUMAN_FEATURE_TO_GENE,
        cols = c("assay", "feature_id"))

# Save
usethis::use_data(HUMAN_FEATURE_TO_GENE, overwrite = TRUE,
                  compress = TRUE, version = 3)
