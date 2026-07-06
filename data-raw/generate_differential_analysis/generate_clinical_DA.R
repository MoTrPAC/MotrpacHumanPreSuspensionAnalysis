# So this function is going to look a bit different than the other generate DA functionality, partially because these clinical chemistry inputs
# were originally categorized as metabolites and some of the samples were processed at different times and such

library(MotrpacHumanPreSuspensionData)
library(MotrpacHumanPreSuspensionAnalysis)
library(dplyr)
library(here)

config_file = "~/config.json"
config = jsonlite::fromJSON(config_file)
repo_local_dir = config$precovid_repo_path
output_local = file.path(repo_local_dir, "data", "tmp", "freeze", "clinical-chemistry", "da")
dir.create(output_local, recursive = TRUE, showWarnings = FALSE)

#for process cov
helpers = file.path(here::here(), "data-raw", "generate_differential_analysis", "generate_differential_modeling_functions.R")
gsutil_helpers = file.path(here::here(), "data-raw", "gsutil_path_parsing.R")
#
source(helpers)
source(gsutil_helpers)

clin_chemistry = MotrpacHumanPreSuspensionData::load_clinical_data()[["chemistry"]]

#note these are all in absolute terms, we log2 transform them and then fit linear models on them pretty directly.

#there is one missing value though. we treat it as missing completely at random.
clin_chem_df = clin_chemistry %>%
  dplyr::bind_rows() %>%
  tibble::column_to_rownames("analyte_name")

#just plug in a random targeted metab platform for metadata and covariates
metab_metadata = MotrpacHumanPreSuspensionData::load_pheno()[["pheno_data"]]

shared_samples = intersect(metab_metadata$vialLabel, colnames(clin_chem_df))

metab_metadata = metab_metadata %>% dplyr::filter(vialLabel %in% shared_samples)

clin_chem_df = clin_chem_df %>%
  dplyr::select(dplyr::any_of(shared_samples)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~log2(.x)))
#make sure they're now matching 1:1 & log transformed, with each row coerced to numeric


process_metadata = process_covariates(meta = metab_metadata,
                                      selected_ome = "metab-t-imm-crt",
                                      tissue_input = "blood",
                                      include_technical = TRUE)

# A few analytes (CK, Cortisol, Glucagon, Insulin) contain scattered missing
# values. The standard run_dream() builds a single contrast matrix L on the full
# sample set, so when dream() drops the NA samples for one of those features the
# design loses columns and L no longer conforms (solve: non-conformable arguments).
# To avoid that, fit each analyte on its own. For every feature we drop the NA
# samples, rebuild the metadata with droplevels() so the design and contrasts
# reflect only the subjects actually measured for that feature, then build a
# feature-specific L and fit. Returns one converted DA table per feature.
fit_clinical_by_feature = function(expression_object,
                                   model_type = "acute",
                                   process_metadata,
                                   tissue,
                                   ome){
  meta_matrix = process_metadata$metadata

  formula = switch(model_type,
                   acute = process_metadata$full_formula,
                   training = process_metadata$training_formula,
                   sex_differences = process_metadata$sex_differences_formula) %>%
    stats::as.formula()

  contrast_generator = switch(model_type,
                              acute = .generate_contrasts_acute,
                              training = .generate_contrasts_training,
                              sex_differences = .generate_sex_contrasts)

  feature_tables = list()
  for(feature in rownames(expression_object)){
    feature_values = unlist(expression_object[feature, ], use.names = TRUE)
    keep_samples = names(feature_values)[!is.na(feature_values)]

    # subset expression and metadata to the subjects actually measured here,
    # then droplevels() so empty group_timepoint / pid levels are removed
    feature_expr = as.matrix(expression_object[feature, keep_samples, drop = FALSE])
    feature_meta = meta_matrix[keep_samples, , drop = FALSE] %>%
      droplevels()
    feature_meta = feature_meta[match(colnames(feature_expr), rownames(feature_meta)), , drop = FALSE]

    # regenerate contrasts on the subset so they reflect the present timepoints/groups
    feature_contrasts = contrast_generator(feature_meta)
    L = variancePartition::makeContrastsDream(formula, feature_meta, contrasts = feature_contrasts)

    fit = variancePartition::dream(feature_expr, formula, feature_meta, L = L)
    fit = variancePartition::eBayes(fit)

    feature_tables[[feature]] = .convert_dream_output(fit,
                                                      metadata = feature_meta,
                                                      tissue = tissue,
                                                      formula = process_metadata$full_formula,
                                                      ome = ome)
  }
  return(dplyr::bind_rows(feature_tables))
}

#only pulling out DA to make sure the columns are in the same order.
da_column_only = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(single_matrix = TRUE)
CLIN_CHEMISTRY_DA = fit_clinical_by_feature(expression_object = clin_chem_df,
                                            model_type = "acute",
                                            process_metadata = process_metadata,
                                            tissue = "blood",
                                            ome = "clinical-chemistry") %>%
  dplyr::mutate(tissue = "blood") %>%
  dplyr::left_join(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER, by = "contrast") %>%
  dplyr::select(dplyr::any_of(colnames(da_column_only)))

file_type = ".txt"
all_file_header = "human-precovid-sed-adu" #this is the base structure for all files within the phase.
tissue_code = "t02-plasma"

#DA table
file_name = paste(all_file_header, tissue_code, "clinical-chemistry", "da", "dream-acute", sep = "_")
file_name = paste0(output_local, "/", file_name, "_v", "1.4", file_type)
write.table(CLIN_CHEMISTRY_DA, file = file_name, row.names = FALSE, sep = '\t', quote = F)


usethis::use_data(CLIN_CHEMISTRY_DA,
                  overwrite = TRUE,
                  version = 3)




