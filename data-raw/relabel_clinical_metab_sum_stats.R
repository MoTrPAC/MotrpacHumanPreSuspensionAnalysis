#!/usr/bin/env Rscript
# One-off relabel, run once for v2.0.1: BLOOD_METAB_T_CLINICAL_SUM_STATS onto the
# metabolomics labelling every other tier already uses — assay = "metab" with the platform in
# its own column.
#
# NOT a generator. Like stack_metab_sum_stats.R beside it, this reads the bytes the package
# already shipped, moves a value between two columns and writes them back. No statistic is
# recomputed. It is idempotent: an object already carrying assay = "metab" is left alone.
#
# Why this way round. The freeze labels every metabolomics DA table assay = "metab" and names
# the platform separately, including the clinical one; BLOOD_METAB_T_CLINICAL_DA ships that
# pair, and precovid-repro's da_assemble_tests.R asserts it by name. The summary-statistic
# side was the odd one out, naming the platform in `assay`, and load_differential_analysis()
# carried a read-time relabel of the DA to paper over the disagreement. Moving this object
# onto the DA's convention is what lets that relabel be deleted.
#
# Being its own object and being labelled "metab" are not in tension. metab-t-clinical stays
# out of the per-tissue stack because it shares five analytes with the research platforms
# (Cortisol, Glycerol, KET, NEFA, Glucose) that would otherwise sit in one object twice;
# `platform` is what tells the two apart, here exactly as on the DA tier.

OBJECT <- "BLOOD_METAB_T_CLINICAL_SUM_STATS"
PLATFORM <- "metab-t-clinical"

# The metabolomics schema, as the stacked objects and the *_DA objects order it.
SCHEMA <- c("tissue", "assay", "platform", "randomGroupCode", "Timepoint",
            "feature_id", "Count", "Mean", "SD")

data_dir <- file.path(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")), "data")
path <- file.path(data_dir, paste0(OBJECT, ".rda"))

e <- new.env()
load(path, envir = e)
x <- as.data.frame(e[[OBJECT]])

if (identical(unique(as.character(x$assay)), "metab")) {
  message(OBJECT, " already carries assay = \"metab\" — nothing to do")
} else {
  stopifnot(identical(unique(as.character(x$assay)), PLATFORM))

  x$platform <- PLATFORM
  x$assay <- "metab"
  x <- x[, SCHEMA, drop = FALSE]
  rownames(x) <- NULL

  assign(OBJECT, x)
  save(list = OBJECT, file = path, version = 3)
  tools::resaveRdaFiles(path, compress = "auto")

  message(sprintf("%s: %d rows, assay = \"%s\", platform = \"%s\"",
                  OBJECT, nrow(x), unique(x$assay), unique(x$platform)))
}
