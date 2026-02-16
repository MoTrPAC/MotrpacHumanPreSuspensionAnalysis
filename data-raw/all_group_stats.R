#We aren't able to offer the full sample level expression data as publically available,
#so instead I'm just generating sum stats for group/timepoint combinations
library(MotrpacHumanPreSuspensionData)
library(dplyr)
library(tidyr)
library(tibble)
library(here)

repo_local_dir = "~/Downloads"
qc_data = MotrpacHumanPreSuspensionData::load_qc(selected_omes = "epigen-atac-seq",
                                                 epigen = TRUE,
                                                 repo_local_dir = repo_local_dir)

da_data = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(epigen = TRUE,
                                                                        single_matrix = TRUE)

# qc_data = MotrpacHumanPreSuspensionData::load_qc(epigen = FALSE,
#                                                  repo_local_dir = repo_local_dir)

.save_one = function(obj, name) {
  assign(x = name, value = obj)

  do.call(
    usethis::use_data,
    args = list(as.name(name), overwrite = TRUE, version = 3)
  )
}

#older code. not really efficient but should find sum stats easily.

for(tissue in names(qc_data)){
  for(assay in names(qc_data[[tissue]])){
    curr_data = qc_data[[tissue]][[assay]][["qc_norm"]]
    if(nrow(curr_data) == 0) next
    curr_meta = qc_data[[tissue]][[assay]][["sample_metadata"]] %>%
      filter(visitcode == "ADU_BAS")
    curr_data = curr_data[,colnames(curr_data) %in% curr_meta$vialLabel]

    curr_meta = curr_meta[match(colnames(curr_data), curr_meta$vialLabel), ] #reorder so they're in the same order as colnames of raw_counts
    vial_to_participant = setNames(curr_meta$pid, curr_meta$vialLabel)
    vial_to_group <- setNames(curr_meta$randomGroupCode, curr_meta$vialLabel)
    vial_to_timepoint <- setNames(curr_meta$Timepoint, curr_meta$vialLabel)

    sample_participant = vial_to_participant[colnames(curr_data)]
    sample_groups <- vial_to_group[colnames(curr_data)]
    sample_timepoints <- vial_to_timepoint[colnames(curr_data)]

    curr_data_long <- as.data.frame(t(curr_data))
    curr_data_long$Sample <- rownames(curr_data_long)
    curr_data_long <- pivot_longer(curr_data_long, -Sample, names_to = "feature_id", values_to = "Value") %>%
      filter(!is.na(Value))
    curr_data_long$randomGroupCode <- sample_groups[curr_data_long$Sample]
    curr_data_long$Timepoint <- sample_timepoints[curr_data_long$Sample]
    curr_data_long$Participant = sample_participant[curr_data_long$Sample]

    group_stats <- curr_data_long %>%
      group_by(randomGroupCode, feature_id, Timepoint) %>%
      summarize(Count = n(),
                Mean = mean(Value),
                SD = sd(Value)
      )

    group_stats$tissue = tissue
    group_stats$assay = assay

    object_name = paste(toupper(tissue), toupper(assay), "SUM_STATS", sep = "_")
    object_name = gsub("-", "_", object_name)
    .save_one(group_stats, name = object_name)
  }
}
