#' Generate QC-normalized metabolomics datasets for Pre-COVID analysis
#'
#' Loads raw metabolomics data from local storage, removes predefined outliers,
#' re-runs KNN imputation, log2 transformation, and median-MAD normalization
#' (untargeted only), selects the appropriate normalized layer per dataset, and
#' writes qc_norm, sample_metadata, and feature_metadata files using
#' \code{write_with_path_name}.
#'
#' This function consolidates the logic from \code{generate_metab_qc_part2.Rmd}
#' into a callable form matching the pattern of other \code{generate_*_qc_norm}
#' functions. Raw data must already be present in
#' \code{file.path(repo_local_dir, "data", "tmp", "metabolomics_qc_norm")}.
#'
#' A JSON config file is required because this function sources tools from two
#' external repositories (\code{motrpac_bic_norm_qc} and the pre-COVID analyses
#' repo). The config must contain at minimum \code{motrpac_bic_norm_qc_repo_path}
#' and \code{precovid_repo_path}. See \code{example_config.json} in the repo root
#' for the expected structure.
#'
#' @param config_file Character scalar. Path to the JSON config file. Defaults
#'   to \code{"~/config.json"}.
#'
#' @return Called for its side effects. Writes files under
#'   \code{file.path(precovid_repo_path, "data", "tmp", "metab_outliers")}.
#'   Invisibly returns the final \code{metabolomics_processed_datasets} list.
#'
#' @keywords internal metabolomics normalization
#' @author christopher jin

generate_metab_qc_norm = function(config_file = "~/config.json") {
  if(!file.exists(config_file))
    stop(config_file, " not found. This function requires a config file to locate motrpac_bic_norm_qc_repo_path and precovid_repo_path.")

  config = jsonlite::fromJSON(config_file)
  repo_local_dir = config$precovid_repo_path

  source(file.path(normalizePath(config$motrpac_bic_norm_qc_repo_path), "tools/unsupervised_normalization_functions.R"))
  source(file.path(normalizePath(config$motrpac_bic_norm_qc_repo_path), "tools/MetabolomicsLibrary.R"))
  source(file.path(normalizePath(config$motrpac_bic_norm_qc_repo_path), "tools/metabolomics_data_parsing_functions.R"))
  source(file.path(normalizePath(config$motrpac_bic_norm_qc_repo_path), "tools/gcp_functions.R"))
  source(file.path(normalizePath(config$precovid_repo_path), "library/JF-shared-functions.R"))
  source(file.path(here(), "data-raw", "gsutil_path_parsing.R"))

  imputation_alpha = 0.2
  max_allowed_NA_rate = 0.2
  required_meta_vars = c("codedsiteid", "visitcode", "pid", "BID", "Timepoint", "randomGroupCode",
                          "htcmavg_hwwt", "wtkg_pcaa", "calculatedAge", "Sex")

  # Load raw metabolomics data
  metab_files = file.path(repo_local_dir, "data", "tmp", "metabolomics_qc_norm")
  metab_obj = list(
    downloaded_files = list.files(metab_files, full.names = TRUE, recursive = TRUE),
    local_path = metab_files
  )
  metab_datasets = suppressWarnings(read_metabolomics_datasets_from_download_obj(metab_obj, TRUE))
  metabolomics_parsed_datasets = metab_datasets$metabolomics_parsed_datasets

  for(currname in c("t02-plasma,metab-t-conv,named", "t02-plasma,metab-t-imm-crt,named")) {
    metabolomics_parsed_datasets[[currname]]$sample_meta$site = "duke"
  }

  pheno_data = MotrpacHumanPreSuspensionData::load_pheno(load_acute_only = FALSE)
  outliers = MotrpacHumanPreSuspensionAnalysis::OUTLIERS
  outliers$vialLabel = as.character(outliers$vialLabel)

  # Filter loop: apply outlier removal and sample filtering per parsed dataset
  for(currname in names(metabolomics_parsed_datasets)) {
    currname_for_outlier = gsub(",unnamed", "", currname)
    currname_for_outlier = gsub(",named", "", currname_for_outlier)
    currname_for_outlier = strsplit(currname_for_outlier, ",")[[1]][2]

    specific_outliers = outliers %>%
      dplyr::filter(stringr::str_starts(ome, currname_for_outlier)) %>%
      dplyr::pull(vialLabel) %>%
      as.character()

    d = metabolomics_parsed_datasets[[currname]]

    d$sample_data = d$sample_data %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("metabolite_name") %>%
      dplyr::select(dplyr::any_of(pheno_data$pheno_data$vialLabel[!is.na(pheno_data$pheno_data$randomGroupCode)])) %>%
      dplyr::select(-dplyr::any_of(specific_outliers))

    if(!is.null(dim(d$control_data)) && ncol(d$control_data) > 1) {
      d$control_data = d$control_data %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames("metabolite_name")
    }

    d$row_annot = .reconstruct_anno(rownames(d$sample_data), d$row_annot)
    d$row_annot$is_named = ifelse(grepl(",named", currname), TRUE, FALSE)
    metabolomics_parsed_datasets[[currname]] = d
  }
  gc()

  # Merge and normalize loop
  metabolomics_processed_datasets = c()
  named_versions = names(metabolomics_parsed_datasets)[-grep("unnamed", names(metabolomics_parsed_datasets))]

  for(currname in named_versions) {
    unnamed_version = gsub("named", "unnamed", currname)

    if(unnamed_version %in% names(metabolomics_parsed_datasets)) {
      new_processed_dataset = merge_named_and_unnamed_metabolomics_datasets(
        metabolomics_parsed_datasets[[currname]],
        metabolomics_parsed_datasets[[unnamed_version]],
        strict = FALSE
      )
    } else {
      new_processed_dataset = metabolomics_parsed_datasets[[currname]]
    }

    new_processed_dataset$sample_meta$BID = substr(new_processed_dataset$sample_meta$sample_id, start = 1, stop = 5)

    curr_meta = pheno_data$pheno_data %>%
      dplyr::filter(vialLabel %in% colnames(new_processed_dataset$sample_data))
    curr_meta$sample_order = new_processed_dataset$sample_meta[rownames(curr_meta), "sample_order"]
    curr_meta$sample_type = new_processed_dataset$sample_meta[rownames(curr_meta), "sample_type"]
    curr_meta$raw_intensity_pre = log2(apply(new_processed_dataset$sample_data[, rownames(curr_meta)], 2, mean, na.rm = TRUE) + 1)

    swap_names = .organize_some_names(currname, new_processed_dataset, metabolomics_processed_datasets)
    newname = swap_names[["newname"]]
    metabolomics_processed_datasets = swap_names[["metabolomics_processed_datasets"]]

    curr_data = new_processed_dataset$sample_data
    if(any(curr_data == 0 | curr_data < 0, na.rm = TRUE)) {
      curr_data[curr_data == 0 | curr_data < 0] = NA
    }
    curr_meta$num_NAs = colSums(is.na(curr_data[, rownames(curr_meta)]))

    impute_as_needed = .organize_missingness(curr_data, new_processed_dataset,
                                              metabolomics_processed_datasets, newname,
                                              max_allowed_NA_rate, imputation_alpha)
    curr_data_imp = impute_as_needed[["curr_data_imp"]]
    metabolomics_processed_datasets = impute_as_needed[["metabolomics_processed_datasets"]]

    log2_imputed = log2(1 + curr_data_imp)
    curr_meta$raw_intensity_post = log2(apply(curr_data_imp[, rownames(curr_meta)], 2, mean, na.rm = TRUE) + 1)

    metabolomics_processed_datasets[[newname]]$sample_meta = curr_meta
    metabolomics_processed_datasets[[newname]]$normalized_data = list(
      raw_imputed = curr_data_imp,
      raw = curr_data,
      log2_imputed = log2_imputed
    )

    if(grepl("metab-u", currname)) {
      median_mad_data = run_median_mad_norm(log2_imputed, zero_medians = FALSE)
      metabolomics_processed_datasets[[newname]]$normalized_data$log2_median_mad = median_mad_data
    }

    metabolomics_processed_datasets[[newname]]$row_annot = impute_as_needed[["curr_data_row_annot"]]
    metabolomics_processed_datasets[[newname]]$rows_to_remove = impute_as_needed[["rows_to_remove"]]
    message("Finished: ", newname)
  }

  na_summary = do.call(rbind, lapply(names(metabolomics_processed_datasets), function(dataset_name) {
    raw = metabolomics_processed_datasets[[dataset_name]][["normalized_data"]][["raw"]]
    #all possible rows to remove -> then we actually check if the removal is "true"
    rows_to_remove = metabolomics_processed_datasets[[dataset_name]][["rows_to_remove"]]
    features_to_exclude = if(is.logical(rows_to_remove)) names(rows_to_remove)[rows_to_remove] else character(0)
    tissue = strsplit(dataset_name, ",")[[1]][1]
    platform = strsplit(dataset_name, ",")[[1]][2]

    data.frame(
      tissue = tissue,
      platform = platform,
      feature_id = rownames(raw),
      pct_na_imputed = rowMeans(is.na(raw))*100
    ) %>%
      dplyr::filter(!feature_id %in% features_to_exclude)
  }))


  # KW test loop: determine sample-centering decision for untargeted datasets
  kw_intensity_report = c()
  for(currname in names(metabolomics_processed_datasets)) {
    if(!grepl("metab-u", currname)) next

    curr_data = metabolomics_processed_datasets[[currname]]
    test_ss = data.frame(
      vialLabel = colnames(curr_data$normalized_data$raw_imputed),
      samp.med = apply(curr_data$normalized_data$raw_imputed, 2, median, na.rm = TRUE),
      quant.75 = apply(curr_data$normalized_data$raw_imputed, 2, function(x) quantile(x, 0.75, na.rm = TRUE))
    )

    curr_int_data = metabolomics_processed_datasets[[currname]]$sample_meta %>%
      dplyr::select(dplyr::any_of(required_meta_vars)) %>%
      tibble::rownames_to_column("vialLabel") %>%
      dplyr::right_join(., test_ss, by = "vialLabel")

    if(length(unique(curr_int_data$Sex)) > 1) {
      sex_pval_1 = kruskal.test(samp.med ~ Sex, data = curr_int_data)$p.value
      sex_pval_2 = kruskal.test(quant.75 ~ Sex, data = curr_int_data)$p.value
    } else {
      sex_pval_1 = NA_real_
      sex_pval_2 = NA_real_
    }

    sex_intensity_report_temp = data.frame(
      dataset = currname,
      group = "Sex",
      outcome = c("samp.med", "quant.75"),
      p.value = c(sex_pval_1, sex_pval_2)
    )

    group_intensity_report_temp = c()
    for(s in unique(curr_int_data$Sex)) {
      curr_s_data = curr_int_data %>% dplyr::filter(Sex == s)
      group_pval_1 = kruskal.test(samp.med ~ randomGroupCode, data = curr_s_data)$p.value
      group_pval_2 = kruskal.test(quant.75 ~ randomGroupCode, data = curr_s_data)$p.value
      group_intensity_report_temp = rbind(group_intensity_report_temp,
                                          data.frame(
                                            dataset = currname,
                                            group = sprintf("randomGroupCode,%s", s),
                                            outcome = c("samp.med", "quant.75"),
                                            p.value = c(group_pval_1, group_pval_2)
                                          ))
    }
    kw_intensity_report = rbind(kw_intensity_report, sex_intensity_report_temp, group_intensity_report_temp)
  }

  sample_centering_decision = kw_intensity_report %>%
    dplyr::group_by(dataset) %>%
    dplyr::slice_min(p.value) %>%
    dplyr::mutate(sample_ctr = ifelse(p.value > 0.01, 1, 0)) %>%
    dplyr::slice_min(sample_ctr) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::filter(sample_ctr == 1) %>%
    dplyr::pull(dataset)

  # Data selection loop: assign data_use per dataset
  for(currname in names(metabolomics_processed_datasets)) {
    curr_norm_data = metabolomics_processed_datasets[[currname]]$normalized_data

    if(grepl("metab-t", currname)) {
      if(nrow(curr_norm_data$raw) > 12) {
        curr_data_use = curr_norm_data$log2_imputed
      } else {
        curr_data_use = curr_norm_data$raw
        for(i in 1:nrow(curr_data_use)) {
          if(sum(is.na(curr_data_use[i,])) != 0) {
            curr_data_use[is.na(curr_data_use)] = min(curr_data_use[i,], na.rm = TRUE) / 2
          }
        }
        curr_data_use = log2(1 + curr_data_use)
      }
    } else if(grepl("metab-u", currname)) {
      if(currname %in% sample_centering_decision) {
        curr_data_use = curr_norm_data$log2_median_mad
      } else {
        curr_data_use = curr_norm_data$log2_imputed
      }
    } else {
      message("Dataset not targeted or untargeted, skipping: ", currname)
      next
    }
    metabolomics_processed_datasets[[currname]]$data_use = curr_data_use
  }

  # Write loop
  outdir = file.path(repo_local_dir, "data", "tmp", "freeze")
  for(dataset in names(metabolomics_processed_datasets)) {
    if(is.null(metabolomics_processed_datasets[[dataset]]$data_use)) next

    ome = strsplit(dataset, ",")[[1]][2]

    tissue_code = strsplit(dataset, ",")[[1]][1]
    tissue = OME_TISSUE_CODE %>% filter(tissue_code == !!tissue_code, ome == !!ome) %>% pull(tissue)

    if(grepl("metab-t", ome)) {
      single_outdir = file.path(outdir, "metabolomics-targeted")
    } else {
      single_outdir = file.path(outdir, "metabolomics-untargeted")
    }
    dir.create(file.path(single_outdir, "qc-norm"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(single_outdir, "metadata"), recursive = TRUE, showWarnings = FALSE)

    qc_norm_table = metabolomics_processed_datasets[[dataset]]$data_use %>%
      as.data.frame() %>%
      dplyr::mutate(feature_id = rownames(.)) %>%
      dplyr::select(feature_id, dplyr::everything()) %>%
      dplyr::filter(!stringr::str_detect(rownames(.), "(?i)istd|standard"))

    metadata_samples = metabolomics_processed_datasets[[dataset]]$sample_meta

    relevant_impute_table = na_summary %>%
      dplyr::filter(tissue == tissue_code, platform == ome) %>%
      dplyr::select(feature_id, pct_na_imputed)

    metadata_features = metabolomics_processed_datasets[[dataset]]$row_annot %>%
      as.data.frame() %>%
      dplyr::rename(feature_id = metabolite_name) %>%
      dplyr::filter(feature_id %in% qc_norm_table$feature_id) %>%
      dplyr::left_join(., relevant_impute_table, by = "feature_id")

    write_with_path_name(qc_norm_table,
                         local_path = file.path(single_outdir, "qc-norm/"),
                         tissue = tissue,
                         ome = ome,
                         data_category = "qc-norm",
                         data_details = "log2",
                         version = "1.4")

    write_with_path_name(metadata_samples,
                         local_path = file.path(single_outdir, "metadata/"),
                         tissue = tissue,
                         ome = ome,
                         data_category = "metadata",
                         data_details = "samples",
                         version = "1.4")

    write_with_path_name(metadata_features,
                         local_path = file.path(single_outdir, "metadata/"),
                         tissue = tissue,
                         ome = ome,
                         data_category = "metadata",
                         data_details = "features",
                         version = "1.4")

    message("Written: ", dataset)
  }

  invisible(metabolomics_processed_datasets)
}



# Private: reorders annotation to match sample data row order
.reconstruct_anno = function(metabolites, anno) {
  if(all(metabolites %in% rownames(anno))) return(anno[metabolites,])
  m = c()
  for(metab in metabolites) {
    ind = which(anno[,1] == metab)
    if(length(ind) < 1) {
      print(paste0("ERROR in .reconstruct_anno: metabolite ", metab, " does not appear in the first column"))
      return(NULL)
    }
    m = rbind(m, anno[ind[1],])
    rownames(m)[nrow(m)] = metab
  }
  return(m)
}

# Private: parses tissue/platform/site and initializes the processed dataset list entry
.organize_some_names = function(currname, new_processed_dataset, metabolomics_processed_datasets) {
  curr_tissue = strsplit(currname, split = ",")[[1]][1]
  curr_platform = strsplit(currname, split = ",")[[1]][2]
  curr_site = tolower(unique(new_processed_dataset$sample_meta$site))
  split_site = strsplit(curr_site, split = " ")[[1]]
  short_site_name = ifelse(length(split_site) > 0,
                           ifelse(split_site[1] == "university", split_site[length(split_site)], split_site[1]),
                           "")
  newname = paste(curr_tissue, curr_platform, short_site_name, sep = ",")

  metabolomics_processed_datasets[[newname]] = list()
  metabolomics_processed_datasets[[newname]]$tissue = curr_tissue
  metabolomics_processed_datasets[[newname]]$site = curr_site
  metabolomics_processed_datasets[[newname]]$platform = curr_platform
  metabolomics_processed_datasets[[newname]]$is_targeted = grepl("metab-t", curr_platform)
  metabolomics_processed_datasets[[newname]]$control_data = new_processed_dataset$control_data
  metabolomics_processed_datasets[[newname]]$control_sample_type =
    new_processed_dataset$sample_meta[new_processed_dataset$sample_meta$sample_id %in%
                                        colnames(new_processed_dataset$control_data),]
  return(list(newname = newname, metabolomics_processed_datasets = metabolomics_processed_datasets))
}

# Private: feature-level filtering + KNN imputation; returns updated list components
.organize_missingness = function(curr_data, new_processed_dataset,
                                 metabolomics_processed_datasets, newname,
                                 max_allowed_NA_rate, imputation_alpha) {
  if(nrow(curr_data) == 1 || is.null(dim(curr_data))) {
    new_processed_dataset$row_annot$num_NAs =
      rowSums(is.na(curr_data[rownames(new_processed_dataset$row_annot),]))
    return(list(
      curr_data_imp = curr_data,
      metabolomics_processed_datasets = metabolomics_processed_datasets,
      curr_data_row_annot = new_processed_dataset$row_annot,
      rows_to_remove = "One feat - no removing"
    ))
  }

  rows_to_remove = !(apply(curr_data, 1, sd, na.rm = TRUE) > 0) |
    rownames(curr_data) == "" |
    rownames(curr_data) == "-" |
    rownames(curr_data) == "_" |
    is.na(rownames(new_processed_dataset$sample_data))
  rows_to_remove[is.na(rows_to_remove)] = TRUE
  row_percent_na = rowSums(is.na(curr_data)) / ncol(curr_data)
  rows_to_remove = rows_to_remove | row_percent_na > max_allowed_NA_rate

  curr_data = curr_data[!rows_to_remove,]
  curr_data_row_annot = new_processed_dataset$row_annot[!rows_to_remove,]
  curr_data_row_annot$num_NAs = rowSums(is.na(curr_data[rownames(curr_data_row_annot),]))
  metabolomics_processed_datasets[[newname]]$control_data = new_processed_dataset$control_data[!rows_to_remove,]

  if(any(is.na(curr_data))) {
    if(sum(is.na(curr_data)) == 1) {
      na_per_metabolite = setNames(colnames(curr_data)[colSums(is.na(curr_data)) > 0],
                                   rownames(curr_data)[rowSums(is.na(curr_data)) > 0])
      na_per_sample = setNames(rownames(curr_data)[rowSums(is.na(curr_data)) > 0],
                               colnames(curr_data)[colSums(is.na(curr_data)) > 0])
      metabolomics_processed_datasets[[newname]]$na_data =
        list(by_metabolite = na_per_metabolite, by_sample = na_per_sample)
    }
    if(sum(is.na(curr_data)) > 1) {
      na_data_temp = curr_data[!complete.cases(curr_data), colSums(is.na(curr_data)) > 0]
      na_per_metabolite = vector(mode = "list", length = nrow(na_data_temp))
      names(na_per_metabolite) = rownames(na_data_temp)
      for(metabolite in names(na_per_metabolite)) {
        na_per_metabolite[[metabolite]] = colnames(na_data_temp[which(is.na(na_data_temp[metabolite,]))])
      }
      na_per_sample = vector(mode = "list", length = ncol(na_data_temp))
      names(na_per_sample) = colnames(na_data_temp)
      for(labelid in names(na_per_sample)) {
        na_per_sample[[labelid]] = rownames(na_data_temp)[which(is.na(na_data_temp[,labelid]))]
      }
      metabolomics_processed_datasets[[newname]]$na_data =
        list(by_metabolite = na_per_metabolite, by_sample = na_per_sample)
    }
    curr_data_imp = min_val_knn_hybrid_imputation(curr_data, imputation_alpha,
                                                  ml_method = "knn",
                                                  num_cores = min(15, nrow(curr_data)))
  } else {
    curr_data_imp = curr_data
  }

  return(list(
    curr_data_imp = curr_data_imp,
    metabolomics_processed_datasets = metabolomics_processed_datasets,
    curr_data_row_annot = curr_data_row_annot,
    rows_to_remove = rows_to_remove
  ))
}


.fix_refmet_names = function(){

}
