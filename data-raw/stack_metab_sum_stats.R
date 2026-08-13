#!/usr/bin/env Rscript
# One-off reshape, run once for v2.0.1: the per-platform metabolomics *_SUM_STATS objects ->
# one stacked {TISSUE}_METAB_SUM_STATS per tissue.
#
# NOT a generator, and not a second one beside the pipeline — see all_group_stats.R in this
# directory for why the package keeps only one. This recomputes nothing. It reads the bytes
# the package already shipped, moves a column, binds them and deletes its inputs, so it
# cannot run a second time: the objects it reads are gone after the first run and it stops.
# Kept in the tree so the reshape is traceable rather than appearing as 32 deletions and 3
# unexplained additions.
#
# The differential analysis has always been keyed this way — step 10 of the reproduction
# pipeline collapses the research metabolomics platforms into a single {TISSUE}_METAB_DA
# carrying assay = "metab" with the platform in its own column — while the summary statistics
# shipped one object per platform naming the platform in `assay`. A caller loading both got a
# list nested one way on one side and another way on the other, and had to translate between
# the two vocabularies before it could join them.
#
# This rebuilds the shipped objects from the shipped objects: no recomputation, the same rows
# under the layout the DA uses. Going forward the objects come from
# precovid-repro/scripts/10_build_data/11_build_sum_stats/all_group_stats.R, which produces
# exactly this shape.
#
# metab-t-clinical is NOT part of the stack, on either tier. It is clinical chemistry, it is
# published as its own object, and it shares five analytes with the research platforms
# (Cortisol, Glycerol, KET, NEFA, Glucose) that would land in one object twice.

suppressWarnings(suppressMessages(library(dplyr)))

data_dir <- file.path(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")), "data")

# Column order: the *_DA order restricted to the columns the two tiers share — what the rows
# are (tissue, assay, platform), then what identifies a row (group, timepoint, feature), then
# the statistics. DA puts its contrast columns between the two; this tier has none.
SCHEMA <- c("tissue", "assay", "platform", "randomGroupCode", "Timepoint",
            "feature_id", "Count", "Mean", "SD")

read_rda <- function(path) {
  e <- new.env()
  load(path, envir = e)
  as.data.frame(e[[ls(e)[1]]])
}

for (tissue in c("ADIPOSE", "BLOOD", "MUSCLE")) {
  pattern <- sprintf("^%s_METAB_[A-Z_]+_SUM_STATS[.]rda$", tissue)
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  files <- files[basename(files) != "BLOOD_METAB_T_CLINICAL_SUM_STATS.rda"]
  if (!length(files)) stop(tissue, ": no per-platform metabolomics summary statistics found")

  parts <- lapply(sort(files), function(f) {
    object <- sub("[.]rda$", "", basename(f))
    x <- read_rda(f)
    # `assay` named the platform on this tier and only on this tier. It moves to `platform`,
    # and `assay` becomes what the DA calls the stack.
    x$platform <- as.character(x$assay)
    x$assay <- "metab"
    x[, SCHEMA, drop = FALSE]
  })

  stacked <- dplyr::bind_rows(parts)
  rownames(stacked) <- NULL

  # The parts are disjoint in (randomGroupCode, feature_id, Timepoint): the CV collapse the
  # metabolomics feature space goes through keeps each RefMet name on exactly one platform
  # per tissue. Asserted rather than assumed — a collision would make the stack silently
  # ambiguous for every caller keying on those three columns.
  key <- paste(stacked$randomGroupCode, stacked$feature_id, stacked$Timepoint)
  if (anyDuplicated(key))
    stop(tissue, ": ", sum(duplicated(key)),
         " (randomGroupCode, feature_id, Timepoint) key(s) appear on more than one platform")

  object <- paste0(tissue, "_METAB_SUM_STATS")
  out <- file.path(data_dir, paste0(object, ".rda"))
  assign(object, stacked)
  save(list = object, file = out, version = 3)
  # The objects this replaces are bzip2; "auto" picks whichever of the three is smallest for
  # this object rather than assuming that carries over to a table twelve times the size.
  tools::resaveRdaFiles(out, compress = "auto")

  file.remove(files)

  message(sprintf("%-24s %6d rows, %2d platforms, %5d features  (replaced %d object(s))",
                  object, nrow(stacked), dplyr::n_distinct(stacked$platform),
                  dplyr::n_distinct(stacked$feature_id), length(files)))
}
