library(MotrpacHumanPreSuspensionAnalysis)

# Save RefMet metabolite names and use to obtain the chemical subclasses with
# https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmetF_form.php
# (accessed 2024-08-07)
refmet_names <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
  ungroup() %>%
  filter(platform == "metabolomics") %>%
  # Remove internal standards
  filter(!is.na(refmet_name),
         !grepl("[[]iSTD[]]|InternalStandard",
                refmet_name, ignore.case = TRUE)) %>%
  select(refmet_name)

write.table(refmet_names, file = file.path("data-raw", "refmet_names.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

refmet_subclasses <- file.path("data-raw", "refmet_results.txt") %>%
  read.delim() %>%
  distinct(Sub.class, Input.name) %>%
  filter(Sub.class != "") %>%
  {split(x = .[["Input.name"]], f = .[["Sub.class"]])}

names(refmet_subclasses) <- paste0("REFMET_", names(refmet_subclasses))

# Write named list to a GMT file
file <- file.path("data-raw", "gmt_processing", "gmt_files",
                  "metabolomics.workbench.refmet.2024.08.07.metabolites.gmt")

MotrpacHumanPreSuspensionAnalysis:::.writeGMT(x = refmet_subclasses,
                                            path = path)


# Compress file
R.utils::gzip(file)
