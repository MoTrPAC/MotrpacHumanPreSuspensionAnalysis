# So this function is going to look a bit different than the other generate DA functionality, partially because these clinical chemistry inputs
# were originally categorized as metabolites and some of the samples were processed at different times and such

library(MotrpacHumanPreSuspensionData)
library(dplyr)
library(here)

output_local = "~/Downloads/"
#for process cov
helpers = file.path(here(), "data-raw", "generate_differential_analysis", "generate_differential_modeling_functions.R")
gsutil_helpers = file.path(here(), "data-raw", "gsutil_path_parsing.R")

source(helpers)
source(gsutil_helpers)

clin_chemistry = MotrpacHumanPreSuspensionData::load_clinical_data()[["chemistry"]]

#note these are all in absolute terms, we log2 transform them and then fit linear models on them pretty directly.

#there is one missing value though. we treat it as missing completely at random.
log_clin_chem_df = clin_chemistry %>%
  dplyr::bind_rows() %>%
  tibble::column_to_rownames("analyte_name") %>%
  dplyr::mutate(across(everything(), ~log2(.x)))


#just plug in a random targeted metab platform for metadata and covariates
metab_metadata = MotrpacHumanPreSuspensionData::load_pheno()[["pheno_data"]] %>%
  filter(vialLabel %in% colnames(clin_chem_df))

clin_chem_df = log_clin_chem_df %>%
  select(any_of(metab_metadata$vialLabel))
#make sure they're now matching 1:1

process_metadata = process_covariates(meta = metab_metadata,
                                      selected_ome = "metab-t-imm-crt",
                                      tissue_input = "blood",
                                      include_technical = TRUE)

fit = run_dream(expression_object = clin_chem_df,
                acute_vs_training = "acute",
                process_metadata = process_metadata,
                voom = FALSE,
                parallel = FALSE)

#only pulling out DA to make sure the columns are in the same order.
da_column_only = load_differential_analysis(single_matrix = TRUE)
CLIN_CHEMISTRY_DA = .convert_dream_output(fit,
                                          metadata = process_metadata$original_meta,
                                          tissue = "blood",
                                          formula = process_metadata$full_formula,
                                          ome = "clinical_chemistry") %>%
  mutate(tissue = "blood") %>%
  left_join(CONTRAST_CONVERTER, by = "contrast") %>%
  select(any_of(colnames(da_column_only)))

file_type = ".txt"
all_file_header = "human-precovid-sed-adu" #this is the base structure for all files within the phase.
tissue_code = "t02-plasma"

#DA table
file_name = paste(all_file_header, tissue_code, "clinical_chemistry", "da", "dream-acute", sep = "_")
file_name = paste0(output_local, file_name, "_v", "1.3", file_type)
write.table(CLIN_CHEMISTRY_DA, file = file_name, row.names = F, sep = '\t', quote = F)

# log-QC-norm
log_cln_out = log_clin_chem_df %>% tibble::rownames_to_column("feature_id") %>%
  select(feature_id, everything())
file_name = paste(all_file_header, tissue_code, "clinical_chemistry", "qc-norm", "log2-transformed", sep = "_")
file_name = paste0(output_local, file_name, "_v", "1.3", file_type)
write.table(log_cln_out, file = file_name, row.names = F, sep = '\t', quote = F)

#metadata
full_metadata = MotrpacHumanPreSuspensionData::pheno$data %>% filter(vialLabel %in% colnames(log_cln_out))
file_name = paste(all_file_header, tissue_code, "clinical_chemistry", "metadata", "samples", sep = "_")
file_name = paste0(output_local, file_name, "_v", "1.3", file_type)
write.table(full_metadata, file = file_name, row.names = F, sep = '\t', quote = F)


usethis::use_data(CLIN_CHEMISTRY_DA,
                  overwrite = TRUE,
                  version = 3)



