#' Filter omics features by minimum paired sample size across groups and timepoints
#'
#' This function enforces a minimum *paired* sample size criterion for each molecular
#' feature within a given tissue and omics layer. Pairing is defined relative to the
#' presence of a non-missing \code{pre_exercise} measurement for a given
#' participant–feature combination. Measurements at other timepoints are only retained
#' if the same participant has a valid \code{pre_exercise} value for that feature.
#'
#' The filtering is applied *within timepoint*, such that participants lacking a valid
#' pre-exercise measurement are excluded only for the affected feature and timepoint,
#' rather than removing the participant globally. Features are retained if, after
#' enforcing pairing, **every group–timepoint cell** contains at least
#' \code{paired_n} non-missing observations.
#'
#' This function modifies and returns the input \code{qc_data} object by updating:
#' \itemize{
#'   \item \code{qc_norm}: filtered to qualifying features only
#'   \item \code{feature_metadata}: restricted to retained feature IDs
#' }
#'
#' @param qc_data
#' A nested list containing quality-controlled and normalized omics data, structured
#' as \code{qc_data[[tissue]][[ome]]}. Each ome must contain:
#' \itemize{
#'   \item \code{qc_norm}: a feature × sample matrix
#'   \item \code{sample_metadata}: metadata including \code{vialLabel}, \code{pid},
#'         \code{randomGroupCode}, \code{Timepoint}, and \code{visitcode}
#'   \item \code{feature_metadata}: metadata with a feature identifier column \code{id}
#' }
#'
#' @param tissue
#' Character scalar specifying the tissue to filter (must be a valid name in
#' \code{qc_data}).
#'
#' @param ome
#' Character scalar specifying the omics layer (e.g., transcriptomics, proteomics,
#' metabolomics) within the selected tissue.
#'
#' @param paired_n
#' Integer specifying the minimum number of *paired* observations required per
#' group–timepoint cell for a feature to be retained. Default is \code{3}.
#'
#' @return
#' A modified version of \code{qc_data} in which the selected tissue–ome combination
#' contains only features meeting the paired sample size criterion.
#'
#' @details
#' Pairing is evaluated at the participant–feature level using the \code{pre_exercise}
#' timepoint as the reference. Only samples from visit \code{"ADU_BAS"} are considered.
#' Counts are computed after excluding unpaired observations, and features failing
#' the minimum requirement in any group–timepoint combination are removed entirely.
#' In the pre-suspension dataset, only proteomics, phosphoproteomics are affected by this filtering.
#'
#' @assumptions
#' \itemize{
#'   \item \code{pre_exercise} is the designated baseline timepoint
#'   \item Sample identifiers in \code{qc_norm} match \code{vialLabel}
#'   \item Feature IDs in \code{qc_norm} rownames correspond to \code{feature_metadata$id}
#' }
#' @author christopher jin

filter_paired_n = function(qc_data,
                           tissue,
                           ome,
                           paired_n = 3){
  feature_meta = qc_data[[tissue]][[ome]][["feature_metadata"]]
  curr_data = qc_data[[tissue]][[ome]][["qc_norm"]]
  curr_meta = qc_data[[tissue]][[ome]][["sample_metadata"]] %>%
    dplyr::filter(visitcode == "ADU_BAS")

  curr_data = curr_data[,colnames(curr_data) %in% curr_meta$vialLabel]
  curr_meta = curr_meta[match(colnames(curr_data), curr_meta$vialLabel), ] #reorder so they're in the same order as colnames of raw_counts

  vial_to_participant = setNames(curr_meta$pid, curr_meta$vialLabel)
  vial_to_group = setNames(curr_meta$randomGroupCode, curr_meta$vialLabel)
  vial_to_timepoint = setNames(curr_meta$Timepoint, curr_meta$vialLabel)

  sample_participant = vial_to_participant[colnames(curr_data)]
  sample_groups = vial_to_group[colnames(curr_data)]
  sample_timepoints = vial_to_timepoint[colnames(curr_data)]

  #needed to make the matrix long - probably a better way to do it, but should work.
  curr_data_long = as.data.frame(t(curr_data))
  curr_data_long$Sample = rownames(curr_data_long)
  curr_data_long = tidyr::pivot_longer(curr_data_long, -Sample, names_to = "feature_id", values_to = "Value")
  curr_data_long$randomGroupCode = sample_groups[curr_data_long$Sample]
  curr_data_long$Timepoint = sample_timepoints[curr_data_long$Sample]
  curr_data_long$Participant = sample_participant[curr_data_long$Sample]

  #-----here we just remove the unpaired samples because we're not including this in our criteria.
  # we do this at a timepoint by timepoint stage, because we dont want to exclude entire participants for some timepoints.
  # Identify participants with non-NA values at pre-exercise
  has_pre_exercise = curr_data_long %>% filter(Timepoint == "pre_exercise", !is.na(Value))
  # Filter to include:
  # - non-pre-exercise values only if there's a valid pre-exercise for that participant/feature/group
  # since the df is long alrdy, we just remove rows that dont qualify (and then they dont get counted)
  curr_data_long_filtered = curr_data_long %>%
    dplyr::filter((Timepoint != "pre_exercise" & !is.na(Value) &
              paste(feature_id, Participant) %in% paste(has_pre_exercise$feature_id, has_pre_exercise$Participant)))

  #so this part removes rows for other timepoints if there's no pre-ex for this participant.
  group_stats = curr_data_long_filtered %>%
    dplyr::group_by(randomGroupCode, feature_id, Timepoint) %>%
    dplyr::summarize(Count = sum(!is.na(Value)), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = c("randomGroupCode", "Timepoint"),
                values_from = "Count") %>%
    rowwise() %>%
    dplyr::mutate(
      min_count = min(c_across(where(is.numeric)), na.rm = TRUE), #min per group/tp for any cell
      avg_count = mean(c_across(where(is.numeric)), na.rm = TRUE) #avg per group/tp for any cell
    ) %>%
    ungroup()
  #so now the count for non-pre_ex timepoints exclude non-paired.
  qualifying_features = group_stats %>%
    dplyr::filter(min_count >= paired_n) %>%
    pull(feature_id)

  qc_data[[tissue]][[ome]][["qc_norm"]] = curr_data[rownames(curr_data) %in% qualifying_features,]
  qc_data[[tissue]][[ome]][["feature_metadata"]] = feature_meta %>% filter(id %in% qualifying_features)

  return(qc_data)
}
