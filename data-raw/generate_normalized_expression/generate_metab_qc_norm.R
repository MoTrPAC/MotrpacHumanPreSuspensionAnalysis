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

generate_metab_qc_norm = function(config_file = "~/config.json",
                                  redownload = FALSE) {
  if(!file.exists(config_file))
    stop(config_file, " not found. This function requires a config file to locate motrpac_bic_norm_qc_repo_path and precovid_repo_path.")

  config = jsonlite::fromJSON(config_file)
  repo_local_dir = config$precovid_repo_path

  tmp = file.path(normalizePath(config$precovid_repo_path), "data", "tmp", "metabolomics_qc_norm")
  if(redownload){
    message("This function was run with redownload = FALSE. It will fail if you have not downloaded the raw metabolomics files")
    .download_metab_files(tmp = tmp)
  }

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

  # Sanity check: across the named parsed datasets, does any metabolite_name carry
  # more than one distinct refmet_name? Pull metabolite_name -> refmet_name from each
  # named dataset's row_annot and flag any metabolite_name with >1 refmet mapping.
  # named_refmet = names(metabolomics_parsed_datasets) %>%
  #   grep(",named", ., value = TRUE) %>%
  #   purrr::map_dfr(function(currname) {
  #     metabolomics_parsed_datasets[[currname]]$row_annot %>%
  #       as.data.frame() %>%
  #       dplyr::select(dplyr::any_of(c("metabolite_name", "refmet_name"))) %>%
  #       dplyr::mutate(dataset = currname,
  #                     site = paste(sort(unique(metabolomics_parsed_datasets[[currname]]$sample_meta$site)), collapse = ";"))
  #   })
  #
  # metabolite_name_refmet_conflicts = named_refmet %>%
  #   dplyr::distinct(metabolite_name, refmet_name, dataset, site) %>%
  #   dplyr::group_by(metabolite_name) %>%
  #   dplyr::filter(dplyr::n_distinct(refmet_name) > 1) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::arrange(metabolite_name, refmet_name, dataset, site)

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
          na_cols = is.na(curr_data_use[i, ])
          if(any(na_cols)) {
            curr_data_use[i, na_cols] = min(curr_data_use[i, ], na.rm = TRUE) / 2
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

  # Write loop according to bic data structure format. -----------
  outdir = file.path(repo_local_dir, "data", "tmp", "freeze")
  full_refmet_fixed_table = .build_metab_refmet_map(metabolomics_processed_datasets)

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

    relevant_refmet_fixes = full_refmet_fixed_table %>%
      dplyr::filter(dataset == !!dataset) %>%
      dplyr::select(feature_id, refmet_name, refmet_id, kegg_id)

    metadata_features = metabolomics_processed_datasets[[dataset]]$row_annot %>%
      as.data.frame() %>%
      dplyr::filter(!stringr::str_detect(rownames(.), "(?i)istd|standard")) %>%
      dplyr::rename(feature_id = metabolite_name) %>%
      dplyr::filter(feature_id %in% qc_norm_table$feature_id) %>%
      dplyr::select(-refmet_name) %>%
      dplyr::left_join(., relevant_impute_table, by = "feature_id") %>%
      dplyr::left_join(., relevant_refmet_fixes, by = "feature_id") %>%
      dplyr::mutate(assay = "metab") %>%
      dplyr::select(assay, feature_id, refmet_name, refmet_id, kegg_id, everything())


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

#internal function for downloading the metab raw results from the google cloud bucket.
.download_metab_files = function(tmp,
                                 gsutil_cmd = "gsutil"){

  # Bucket structure is: tissue/platform/results/<data files>. For GET, this path should
  # also have a qa_qc directory with the qc_metrics and sample_metadata files.
  targeted_buckets = c("gs://motrpac-data-hub/human-precovid/results/metabolomics-targeted/")
  untargeted_buckets = c("gs://motrpac-data-hub/human-precovid/results/metabolomics-untargeted/")
  clinical_bucket = c(" gs://motrpac-data-hub/human-precovid/phenotype/human-precovid-sed-adu/raw/data_sets/")
  pheno_bucket = "gs://motrpac-data-hub/human-eqc/adult/" #Human eqc data and randomization information

  for(targeted_bucket in targeted_buckets){
    rem_prev = targeted_bucket == targeted_buckets[1] # TRUE for first iteration of for loop, seems unnecessary
    obj = DownloadBucketLocal(targeted_bucket,tmp,gsutil_cmd)
  }
  for(untargeted_bucket in untargeted_buckets){
    obj = DownloadBucketLocal(untargeted_bucket,tmp,gsutil_cmd)
  }
}


# Reorders or rebuilds a row-annotation data frame to match the order of
# `metabolites`. Fast path: if all names are already rownames, subsets directly.
# Slow path: iterates and matches by the first column, taking the first hit when
# duplicates exist. Returns NULL (with a printed error) if any name is absent.
# metabolites: character vector of metabolite names matching rownames of sample data
# anno: data frame with metabolite names as rownames or in the first column
# Returns: anno reordered/subset to align with metabolites, or NULL on failure
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

# Parses tissue, platform, and site from a comma-delimited dataset name and
# initializes the corresponding entry in metabolomics_processed_datasets with
# tissue/site/platform metadata and control data. Site names of the form
# "University of X" are shortened to "X"; single-word names are used as-is.
# currname: comma-delimited string, e.g. "t02-plasma,metab-t-conv,named"
# new_processed_dataset: parsed dataset list (must contain sample_meta$site)
# metabolomics_processed_datasets: accumulator list built across loop iterations
# Returns: list(newname, metabolomics_processed_datasets)
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

# Filters features with zero variance, blank/placeholder rownames, or NA rate
# above max_allowed_NA_rate, records NA patterns by metabolite and sample, then
# runs KNN imputation via min_val_knn_hybrid_imputation. Single-feature datasets
# skip filtering and imputation entirely. Updates metabolomics_processed_datasets
# in place with control data and NA audit info for newname.
# curr_data: feature × sample numeric matrix
# new_processed_dataset: parsed dataset (provides row_annot and sample_data rownames)
# metabolomics_processed_datasets: accumulator list, updated for newname
# newname: key into metabolomics_processed_datasets being populated
# max_allowed_NA_rate: features with row NA fraction above this are dropped
# imputation_alpha: alpha parameter passed to min_val_knn_hybrid_imputation
# Returns: list(curr_data_imp, metabolomics_processed_datasets, curr_data_row_annot, rows_to_remove)
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


# Applies manual corrections to refmet_name to align all values with the current
# RefMet standard. Two correction sources are merged into one pass: feature_id-
# based overrides for cases where the stored refmet_name is wrong/missing, and a
# name-to-name map covering capitalization errors, typos, slash-delimited
# ambiguities, and lab-internal naming conventions. Some labs submitted
# annotations using outdated LIPID MAPS names or lab-internal aliases that
# predate the current RefMet standard.
# metab: data frame with columns feature_id and refmet_name
# Returns: metab with refmet_name corrected in place
.fix_refmet_names = function(metab) {
  # feature_id-based overrides
  metab <- metab %>%
    mutate(
      refmet_name = case_when(
        feature_id == "13-HODE" ~ "13-HODE",
        feature_id == "13 HODE" ~ "13-HODE", #new for 1.4
        feature_id == "FA(20:1)" ~ "Eicosenoic acid", #new for 1.4
        feature_id == "FA(18:2)" ~ "Linoleic acid", #new for 1.4
        feature_id == "FA(22:6)" ~ "Docosahexaenoic acid", #new for 1.4
        TRUE ~ refmet_name
      )
    )

  refmet_name_map <- c(
    "cholesterol sulfate"              = "Cholesterol sulfate",
    "oleamide"                         = "Oleamide",
    "tridecylamine"                    = "Tridecylamine",
    "PGE2 thanolamide"                 = "PGE2 ethanolamide",
    "12(13)-EpMOE"                     = "12(13)-EpOME",
    "Docosatetraenoic aicd"            = "Docosatetraenoic acid",
    "Tetracosenoic aicd"               = "Tetracosenoic acid",
    "Prostagladin"                     = "Prostaglandin", #new for 1.4

    "Chenodeoxycholic acid\\Deoxycholic acid"   = "Deoxycholic acid",
    "Glycocholic acid\\Glycohyocholic acid"     = "Glycocholic acid",

    "N-Lauroylglycine"                 = "NAGly 12:0",
    "N-linoleoylglycine"               = "NAGly 18:2(9Z,12Z)",
    "N-Oleoyl glycine"                 = "NAGly 18:1(9Z)",
    "N-Undecanoylglycine"              = "NAGly 11:0",
    "N-Myristoylglycine"               = "NAGly 14:0",

    "CAR 12:0-OH"                      = "CAR 12:0;OH",
    "CAR 14:0-OH"                      = "CAR 14:0;OH",
    "CAR 14:1-OH"                      = "CAR 14:1;OH",
    "CAR 4:0-OH"                       = "CAR 4:0;OH",

    "DG 16:0_16:0_0:0"                 = "DG 16:0_16:0",
    "DG 16:0_16:1_0:0"                 = "DG 16:0_16:1",
    "DG 16:0_18:1_0:0"                 = "DG 16:0_18:1",
    "DG 18:2_18:2_0:0"                 = "DG 18:2_18:2",
    "DG 18:2_18:3_0:0"                 = "DG 18:2_18:3",

    "13-OxoODE(13-KODE)"               = "13-Oxo-ODE",
    "2-Arachidonoyl Glycerol (2AG)"    = "MG 0:0/20:4/0:0",
    "5,6-DiHET"                        = "5,6-DiHETE",
    "8(9)-DiHET"                       = "8,9-DiHETE",
    "CoA(15:0)_and_CoA(C14:1-OH)"     = "Pentadecanoyl-CoA/Hydroxytetradecenoyl-CoA",
    "CoA(2:0-COOH)_and_CoA(4:0-OH)"   = "Malonyl-CoA/Hydroxybutyryl-CoA",
    "Linoleoyl Ethanolamide (LEA)"     = "Linoleoyl-EA",
    "Oleoyl Ethanolamide (OEA)"        = "Oleoyl-EA",
    "PC(O-33:2)>PC(O-15:0/18:2)"      = "PC O-16:1/20:4",
    "PC(O-36:5)<PC(O-16:1/20:4)"      = "PC O-16:1/20:4",
    "PE(36:4)>(16:0_20:4)"            = "PE 16:0_20:4",
    "PE(38:4)>(PE(18:0_20:4)"         = "PE 18:0_20:4",
    "Stearoyl Ethanolamide (ceramid)"  = "Stearoyl-EA"
  )

  metab <- metab %>%
    mutate(refmet_name = dplyr::recode(refmet_name, !!!refmet_name_map, .default = refmet_name))

  return(metab)
}

# Eric Leslie's implementation of the Metabolomics Workbench RefMet batch API.
# Queries RefMet to standardize metabolite names and retrieve RefMet IDs and
# KEGG IDs. Platform-specific LC suffixes (_hp_, _rp_, _rn_, _in_, _lp_, _ln_)
# are stripped from names before lookup. All names are submitted in a single
# POST request.
# metab: data frame with columns feature_id, refmet_name, dataset
# Returns: data frame with columns feature_id, platform, refmet_name, refmet_id, kegg_id
.annotate_refmet = function(metab){
  tosearch <- "_hp_|_rp_|_rn_|_in_|_lp_|_ln_"
  metab <- metab %>%
    dplyr::mutate(
      lookup_refmet = dplyr::if_else(
        grepl(tosearch, refmet_name),
        gsub("(.*)(_\\w{2}_\\w{1})", "\\1", refmet_name),
        refmet_name
      ),
      lookup_refmet = trimws(lookup_refmet)
    )

  mets <- stringi::stri_join_list(list(metab$lookup_refmet), sep = "\n")
  h <- curl::new_handle()
  curl::handle_setform(h, metabolite_name = mets)
  req <- curl::curl_fetch_memory(
    "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_minID.php",
    handle = h
  )

  refmet_result <- utils::read.table(
    text = rawToChar(req$content),
    header = TRUE,
    na.strings = "-",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    sep = "\t"
  )

  refmet_result[is.na(refmet_result)] <- '-'
  refmet_result[refmet_result == ''] <- '-'

  refmet_annotated <- refmet_result %>%
    dplyr::transmute(
      lookup_refmet = Input.name,
      refmet_name_std = Standardized.name,
      refmet_id = RefMet_ID,
      kegg_id = KEGG_ID
    ) %>%
    dplyr::left_join(
      metab %>% dplyr::select(feature_id, lookup_refmet, dataset),
      by = "lookup_refmet",
      relationship = "many-to-many"
    ) %>%
    dplyr::transmute(
      feature_id,
      platform = "metabolomics",
      refmet_name = dplyr::na_if(refmet_name_std, "-"),
      refmet_id = dplyr::na_if(refmet_id, "-"),
      kegg_id = dplyr::na_if(kegg_id, "-")
    ) %>%
    dplyr::distinct()

  return(refmet_annotated)
}


#using Eric Leslie's implementation of kegg API(s).

# Note Chris: I got rid of eric's method 1 (Load RefMet to KEGG map from Metabolomics Workbench REST API) because we already use metabolomics workbench for the above refmet mapping.
# Calling metab workbench again added 0 new features.
.annotate_kegg_resources = function(final_metab){

  # ==============================================================================
  # APPLY ADDITIONAL KEGG ANNOTATION METHODS
  # ==============================================================================
  cat("\n=== APPLYING ADDITIONAL KEGG ANNOTATION METHODS ===\n")
  cat("Current KEGG coverage:", sum(!is.na(final_metab$kegg_id)), "/", nrow(final_metab),
      "(", round(100 * mean(!is.na(final_metab$kegg_id)), 1), "%)\n\n")

  # Store original KEGG IDs for comparison
  final_metab <- final_metab %>%
    mutate(kegg_id_original = kegg_id)

  # Method 2: Annotate using KEGGREST
  cat("Querying KEGG database via KEGGREST...\n")

  # Get metabolites that still don't have KEGG IDs
  metab_no_kegg = final_metab %>%
    filter(is.na(kegg_id), !is.na(refmet_name)) %>%
    mutate(kegg_id_keggrest = .get_kegg_ids_via_keggrest(refmet_name))

  final_metab_keggrest <- final_metab %>%
    left_join(
      metab_no_kegg %>% select(feature_id, kegg_id_keggrest),
      by = "feature_id",
      relationship = "many-to-many"
    ) %>%
    mutate(kegg_id = if_else(is.na(kegg_id), kegg_id_keggrest, kegg_id)) %>%
    select(-kegg_id_keggrest) %>%
    distinct() %>%
    select(-kegg_id_original)

  final_metab_keggrest = final_metab_keggrest %>%
    dplyr::group_by(feature_id, refmet_id) %>%
    dplyr::filter(!(is.na(kegg_id) & any(!is.na(kegg_id)))) %>%
    dplyr::ungroup()

  cat("Final KEGG coverage:", sum(!is.na(final_metab_keggrest$kegg_id)), "/", nrow(final_metab_keggrest),
      "(", round(100 * mean(!is.na(final_metab_keggrest$kegg_id)), 1), "%)\n\n")

  return(final_metab_keggrest)

}

# Downloads the full KEGG compound list in a single API call via
# KEGGREST::keggList("compound"), parses all synonyms, and resolves
# metabolite_names by exact local match. Avoids per-compound API queries.
# Multiple KEGG IDs for the same name are collapsed with ";".
# metabolite_names: character vector of RefMet-standardized metabolite names
# Returns: character vector aligned to metabolite_names (NA where no match)
.get_kegg_ids_via_keggrest = function(metabolite_names) {
  all_compounds = KEGGREST::keggList("compound")

  compound_df = data.frame(
    kegg_id = stringr::str_remove(names(all_compounds), "cpd:"),
    raw_name = as.character(all_compounds),
    stringsAsFactors = FALSE
  ) %>%
    tidyr::separate_rows(raw_name, sep = ";") %>%
    dplyr::mutate(raw_name = trimws(raw_name)) %>%
    dplyr::group_by(raw_name) %>%
    dplyr::summarise(kegg_id = paste(kegg_id, collapse = ";"), .groups = "drop")

  data.frame(raw_name = metabolite_names, stringsAsFactors = FALSE) %>%
    dplyr::left_join(compound_df, by = "raw_name") %>%
    dplyr::pull(kegg_id)
}


# Assembles the full feature-to-RefMet annotation table across all metabolomics
# platforms. Collects row annotations from each platform's row_annot, applies
# manual RefMet name corrections (.fix_refmet_names), queries the RefMet batch
# API (.annotate_refmet), and fills residual KEGG IDs via KEGGREST
# (.annotate_kegg_resources). Only named features (is_named == TRUE) are
# included; targeted immuno panels lacking is_named are treated as named.
# metabolomics_processed_datasets: named list output from the normalization loop
# Returns: data frame with feature_id, original_refmet_annotation, refmet_name,
#   refmet_id, kegg_id, dataset, and additional row_annot metadata columns
.build_metab_refmet_map = function(metabolomics_processed_datasets){
  metab = bind_rows(lapply(names(metabolomics_processed_datasets), FUN = function(metab_platform) {
    platform = metabolomics_processed_datasets[[metab_platform]][["row_annot"]]
    specific_platform_name = metab_platform

    # Targeted immuno panels don't have is_named; treat as named
    if (!"is_named" %in% names(platform)) platform$is_named <- TRUE

    platform %>%
      dplyr::mutate(
        is_named = as.logical(is_named),
        across(any_of(c("metabolite_name","refmet_name","formula")), as.character),
        dataset = specific_platform_name,
      )
  })) %>%
    dplyr::rename(feature_id = metabolite_name) %>%
    dplyr::filter(is_named == TRUE)

  manual_fixes_done = .fix_refmet_names(metab)
  refmet_annotations = .annotate_refmet(manual_fixes_done) %>%
    dplyr::select(feature_id, refmet_name, refmet_id, kegg_id)

  #------add refmet stuff via API------
  manual_fixes_refmet = manual_fixes_done %>%
    dplyr::rename(original_refmet_annotation = refmet_name) %>%
    dplyr::left_join(., refmet_annotations, by = c("feature_id"), relationship = "many-to-many")

  #------add additional kegg annotations----------
  final_metab_refmet_map = .annotate_kegg_resources(manual_fixes_refmet)

  return(final_metab_refmet_map)
}



