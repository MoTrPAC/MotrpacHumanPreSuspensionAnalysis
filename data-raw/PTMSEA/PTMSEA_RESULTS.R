# Builds PTMSEA_RESULTS from the combined PTM-SEA GCTs in this folder.
# Run from the package root. See README.md for provenance.
library(cmapR)

PTMSEA_RESULTS <- list(
  muscle = parse_gctx("data-raw/PTMSEA/muscle_prot-ph_ptmsea_EE-RE-vs-CON_combined.gct"),
  adipose = parse_gctx("data-raw/PTMSEA/adipose_prot-ph_ptmsea_EE-RE-vs-CON_combined.gct")
)

usethis::use_data(PTMSEA_RESULTS, overwrite = TRUE)
