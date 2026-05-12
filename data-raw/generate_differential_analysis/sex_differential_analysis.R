library(here)
library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacHumanPreSuspensionData)
library(BiocParallel)
library(dplyr)
library(variancePartition)
library(reformulas)

differential_analysis_scripts = list.files(file.path(here(), "data-raw", "generate_differential_analysis"),
                                           full.names = TRUE)
#dont include the clinical DA because that was implemented slightly differently.
differential_analysis_scripts = differential_analysis_scripts[!grepl("clini|sex_", differential_analysis_scripts)]
differential_analysis_scripts = c(differential_analysis_scripts, file.path(here(), "data-raw", "gsutil_path_parsing.R"))


lapply(differential_analysis_scripts, source)

#-----first we find the qualifying omes + tissues to actually do the analysis on------
# repeated code from `MotrpacHumanPreSuspensionAcute` to generate table s1, but split by sex.

#largely speaking, if we want at least n = 3 per sex + grp + tp, we include, for Muscle + Blood:
#  metab (most platforms), transcriptomics

#for adipose, we can choose the 3.5 hour timepoint, but this gets a little more complicated

all_dataset = load_qc(epigen = FALSE,
                      repo_local_dir = repo_local_dir,
                      remove_redundant_metab = TRUE)

split_tissues = unlist(all_dataset, recursive = FALSE)
split_tissue_assay = unlist(split_tissues, recursive = FALSE)
qc_norm_only = split_tissue_assay[grep("qc_norm", names(split_tissue_assay))]

#makes a big list of tissues, assays
samples_per_plat = data.frame(
  list_name = names(qc_norm_only),
  vialLabels = sapply(qc_norm_only, function(x) paste(colnames(x), collapse = ";"))
) %>%
  mutate(
    tissue = sapply(strsplit(list_name, "\\."), `[`, 1),
    assay = sapply(strsplit(list_name, "\\."), `[`, 2)
  ) %>%
  select(-list_name)

#then uses the pheno file to annotate group and timepoint info and split the vector
qc_norm_counts = lapply(seq_len(nrow(samples_per_plat)), function(row) {
  vials_vec = strsplit(samples_per_plat[row, "vialLabels"], ";")[[1]] %>%
    base::trimws()

  pheno$data %>%
    filter(vialLabel %in% vials_vec) %>%
    group_by(randomGroupCode, Timepoint, Sex) %>%
    summarise(
      n = n(),
      vialLabels = paste(vialLabel, collapse = ","),
      .groups = "drop"
    ) %>%
    mutate(
      tissue = samples_per_plat$tissue[row],
      assay = samples_per_plat$assay[row]
    )
}) %>%
  bind_rows() %>%
  select(tissue, assay, randomGroupCode, Sex, Timepoint, n)

qualifying_omes_tissues = qc_norm_counts %>%
  dplyr::filter(randomGroupCode != "ADUControl") %>%
  group_by(tissue, assay) %>%
  mutate(flag_low_n = any(n < 3)) %>%
  filter(!flag_low_n) %>%
  distinct(tissue, assay)

# If we don't include controls - see filter just above. we CAN do adipose.
# If we require >3 for EVERY cell, this limits our conclusions that we can make.

#--------------------------------

for(row in seq(nrow(qualifying_omes_tissues))){
  if(grepl("metab", qualifying_omes_tissues$assay[row])) next
  generate_DA_inputs(repo_local_dir = "~/Downloads/",
                     model_type = "sex_differences",
                     selected_omes = qualifying_omes_tissues$assay[row],
                     selected_tissues = qualifying_omes_tissues$tissue[row],
                     epigen = FALSE,
                     parallel = FALSE)
}

generate_DA_inputs(repo_local_dir = "~/Downloads/",
                   model_type = "sex_differences",
                   selected_omes = "transcript-rna-seq",
                   selected_tissues = "adipose",
                   epigen = FALSE,
                   parallel = FALSE)

