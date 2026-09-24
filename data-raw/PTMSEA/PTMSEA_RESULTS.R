# Builds PTMSEA_RESULTS from the combined PTM-SEA GCTs in this folder, in the
# long shape of CAMERA_RESULTS. Run from the package root. See README.md for
# provenance.
library(cmapR)
library(dplyr)

gct_files <- c(
  muscle = "data-raw/PTMSEA/muscle_prot-ph_ptmsea_EE-RE-vs-CON_combined.gct",
  adipose = "data-raw/PTMSEA/adipose_prot-ph_ptmsea_EE-RE-vs-CON_combined.gct"
)

load("data/CONTRAST_CONVERTER.rda")
load("data/CAMERA_RESULTS.rda")

gct_to_long <- function(tissue, file) {
  gct <- parse_gctx(file)
  rdesc <- gct@rdesc

  lapply(colnames(gct@mat), function(column) {
    # rdesc column names are make.names() of the mat column names
    suffix <- make.names(column)
    overlap <- rdesc[[paste0("Signature.set.overlap.", suffix)]]

    data.frame(
      tissue = tissue,
      contrast_short = sub("^[a-z]+\\.z\\.std_", "", column),
      set = gct@rid,
      set_size = lengths(strsplit(overlap, "|", fixed = TRUE)),
      set_size_DB = as.integer(rdesc$Signature.set.size),
      NES = gct@mat[, column],
      p_value = as.numeric(rdesc[[paste0("pvalue.", suffix)]]),
      adj_p_value = as.numeric(rdesc[[paste0("fdr.pvalue.", suffix)]])
    )
  }) %>%
    bind_rows()
}

PTMSEA_RESULTS <- Map(gct_to_long, names(gct_files), gct_files) %>%
  bind_rows() %>%
  left_join(select(CONTRAST_CONVERTER, contrast, contrast_short, contrast_type),
            by = "contrast_short") %>%
  mutate(
    tissue = factor(tissue, levels = levels(CAMERA_RESULTS$tissue)),
    assay = factor("prot-ph", levels = levels(CAMERA_RESULTS$assay)),
    contrast_short = factor(contrast_short,
                            levels = levels(CONTRAST_CONVERTER$contrast_short)),
    collection = factor("PTMSIGDB"),
    database = factor(sub("_.*", "", set)),
    set_id = factor(set),
    set_short = paste0(gsub("_", " ", sub("^[^_]+_", "", set)),
                       " (", database, ")"),
    set = factor(set),
    set_short = factor(set_short),
    size_ratio = set_size / set_size_DB,
    direction = factor(ifelse(NES > 0, "Up", "Down"), levels = c("Up", "Down"))
  ) %>%
  select(tissue, assay, contrast_type, contrast, contrast_short, collection,
         database, set_id, set, set_short, set_size, set_size_DB, size_ratio,
         direction, NES, p_value, adj_p_value) %>%
  arrange(tissue, contrast, p_value)

stopifnot(
  !anyNA(PTMSEA_RESULTS$contrast),
  !anyNA(PTMSEA_RESULTS$NES),
  nrow(PTMSEA_RESULTS) == 506 * 6 + 437 * 2,
  !anyDuplicated(PTMSEA_RESULTS[, c("tissue", "contrast", "set_id")]),
  !anyDuplicated(unique(PTMSEA_RESULTS[, c("set_id", "set_short")])$set_short)
)

usethis::use_data(PTMSEA_RESULTS, overwrite = TRUE)
