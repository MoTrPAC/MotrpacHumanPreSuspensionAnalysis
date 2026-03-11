library(here)
library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacHumanPreSuspensionData)
library(BiocParallel)

differential_analysis_scripts = list.files(file.path(here(), "data-raw", "generate_differential_analysis"),
                                           full.names = TRUE, include.dirs = FALSE, no.. = TRUE)
differential_analysis_scripts = differential_analysis_scripts[!grepl("clini|sex_", differential_analysis_scripts)]
differential_analysis_scripts = c(differential_analysis_scripts, file.path(here(), "data-raw", "gsutil_path_parsing.R"))


lapply(differential_analysis_scripts, source)

generate_DA_inputs(repo_local_dir = "~/Downloads/",
                   model_type = "sex_differences",
                   selected_omes = metab_only_list(),
                   epigen = FALSE,
                   parallel = TRUE)


#
# test123 = read.csv("/Users/ch57584/Downloads/data/tmp/freeze_DA/human-precovid-sed-adu_t06-muscle_metab-u-rpneg_da_dream-sex_differences_v1.2.txt",
#                    sep = "\t")

# all_dataset = load_qc(epigen = FALSE, #toggle true if needed
#                       repo_local_dir = "~/Downloads/",
#                       remove_redundant_metab = TRUE)
#
# split_tissues = unlist(all_dataset, recursive = FALSE)
# split_tissue_assay = unlist(split_tissues, recursive = FALSE)
# qc_norm_only = split_tissue_assay[grep("qc_norm", names(split_tissue_assay))]
#
# #makes a big list of tissues, assays
# samples_per_plat = data.frame(
#   list_name = names(qc_norm_only),
#   vialLabels = sapply(qc_norm_only, function(x) paste(colnames(x), collapse = ";"))
# ) %>%
#   mutate(
#     tissue = sapply(strsplit(list_name, "\\."), `[`, 1),
#     assay = sapply(strsplit(list_name, "\\."), `[`, 2)
#   ) %>%
#   select(-list_name)
#
# #then uses the pheno file to annotate group and timepoint info and split the vector
# qc_norm_counts = lapply(seq_len(nrow(samples_per_plat)), function(row) {
#   vials_vec = strsplit(samples_per_plat[row, "vialLabels"], ";")[[1]] %>%
#     base::trimws()
#
#   pheno$data %>%
#     filter(vialLabel %in% vials_vec) %>%
#     group_by(randomGroupCode, Timepoint) %>%
#     summarise(
#       n = n(),
#       vialLabels = paste(vialLabel, collapse = ","),
#       .groups = "drop"
#     ) %>%
#     mutate(
#       tissue = samples_per_plat$tissue[row],
#       assay = samples_per_plat$assay[row]
#     )
# }) %>%
#   bind_rows() %>%
#   select(tissue, assay, randomGroupCode, Timepoint, n, vialLabels)
