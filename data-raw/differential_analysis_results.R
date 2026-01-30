library(dplyr)
library(data.table)

# Extra information about each contrast
CONTRAST_CONVERTER <- file.path("data-raw", "contrast_converter.txt") %>%
  read.delim() %>%
  mutate(
    contrast_type = case_when(
      grepl("^[^.]+\\.pre_exercise - [^.]+\\.pre_exercise$",
            contrast_short) ~ "baseline",
      grepl("^Control", contrast_short) ~ "control_only",
      grepl("Endur", contrast_short) &
        grepl("Resist", contrast_short) ~ "Endur_vs_Resist",
      grepl("^[^. ]+.[^ ]+ - [^.]+\\.[^ ]+$",
            contrast_short) ~ "exercise_no_controls",
      TRUE ~ "exercise_with_controls"
    ),
    contrast_type = factor(contrast_type,
                           levels = c("exercise_with_controls",
                                      "exercise_no_controls",
                                      "Endur_vs_Resist",
                                      "baseline",
                                      "control_only")),
    contrast_order = 1:n(),
    contrast_category = sub("([^.]+)\\.[^-]+ - ([^.]+).*",
                            "\\1-\\2",
                            contrast_short),
    contrast_category = gsub("Control", "CON", contrast_category),
    contrast_category = gsub("Endur", "EE", contrast_category),
    contrast_category = gsub("Resist", "RE", contrast_category),
    across(.cols = c(contrast, contrast_short,
                     contrast_type, contrast_category),
           .fns = ~ factor(.x, levels = unique(.x)))
  ) %>%
  dplyr::mutate(contrast_left = stringr::str_split_fixed(contrast_short, " - ", 2)[,1],
                randomGroupCode = case_when(
                  contrast_category == "Endur_vs_Resist" ~ "ADUEndur - ADUResist",
                  TRUE ~ paste0("ADU", stringr::str_split_fixed(contrast_left, "\\.", 2)[,1])
                ),
                Timepoint = as.factor(stringr::str_split_fixed(contrast_left, "\\.", 2)[,2])) %>%
  dplyr::select(-contrast_left) %>%
  relocate(contrast_order, .before = everything()) %>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("pre_exercise",
                                       "during_20_min",
                                       "during_40_min",
                                       "post_10_min",
                                       "post_15_30_45_min",
                                       "post_3.5_4_hr",
                                       "post_24_hr")))
  #make sure the timepoints are in order of actual timepoints, instead of otherwise

setDT(CONTRAST_CONVERTER)

# Save CONTRAST_CONVERTER
usethis::use_data(CONTRAST_CONVERTER, overwrite = TRUE, version = 3)


## Prepare differential analysis results ----
repo_local_dir <- tempdir()

# Very slow
ls <- MotrpacHumanPreSuspensionAnalysis:::.load_differential_analysis(
  repo_local_dir = repo_local_dir,
  gsutil = "gsutil", # gsutil was added to PATH
  epigen = FALSE
)

# Stack metabolomics platforms by tissue
ls2 <- lapply(ls, function(tissue_i) {
  tissue_i["metab-meta-reg"] <- NULL # remove meta-regression results

  is_metab <- grep("metab", names(tissue_i))

  if (length(is_metab)) {
    metab_df <- bind_rows(tissue_i[is_metab], .id = "platform") %>%
      mutate(assay = "metabolomics") %>%
      relocate(platform, .after = assay)

    tissue_i[is_metab] <- NULL
    tissue_i["metab"] <- list(metab_df)
  }

  return(tissue_i)
})

# Combine tissues and omes
DA_list <- unlist(ls2, recursive = FALSE) %>%
  MotrpacHumanPreSuspensionAnalysis:::.process_raw_DA()

## Rename DA_list elements to create .rda files
tissues <- sub("\\..*", "", names(DA_list))
tissues <- factor(tissues,
                  levels = c("adipose", "blood", "muscle"),
                  labels = c("ADIPOSE", "BLOOD", "MUSCLE"))
tissues <- as.character(tissues)

omes <- sub(".*\\.", "", names(DA_list))
omes <- factor(omes,
               levels = c("transcript-rna-seq",
                          "prot-pr", "prot-ph",
                          "prot-ol", "metab"),
               labels = c("TRNSCRPT", "PROT_PR", "PROT_PH",
                          "PROT_OL", "METAB"))
omes <- as.character(omes)

names(DA_list) <- paste0(tissues, "_", omes, "_DA")


## Save individual .rda objects ---
purrr::walk2(.x = DA_list, .y = names(DA_list), .f = function(obj, name) {
  assign(x = name, value = obj)

  do.call(usethis::use_data,
          args = list(as.name(name), overwrite = TRUE, version = 3))
})


# Remove temporary directory
unlink(repo_local_dir, recursive = TRUE, force = TRUE)
