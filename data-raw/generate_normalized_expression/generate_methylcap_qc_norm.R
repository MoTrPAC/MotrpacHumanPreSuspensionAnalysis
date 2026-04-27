#' Generate feature metadata for MethylCap-seq
#'
#' Extracts feature metadata from the existing \code{MotrpacHumanPreSuspensionData}
#' package object for the \code{epigen-methylcap-seq} assay. Raw MethylCap-seq
#' data are **not processed in R** — all preprocessing, normalization, and
#' quality control are performed using an external, assay-specific pipeline.
#' Only the feature metadata file is generated here.
#'
#' @param config_file Character scalar. Path to the JSON config file. Defaults
#'   to \code{"~/config.json"}.
#'
#' @return Called for its side effects. Writes feature metadata files derived
#'   from the existing package data.

generate_methylcap_qc_norm = function(){
  if(!file.exists(config_file))
    stop(config_file, " not found. This function requires a config file to locate motrpac_bic_norm_qc_repo_path and precovid_repo_path.")

  config = jsonlite::fromJSON(config_file)
  repo_local_dir = config$precovid_repo_path

  data = MotrpacHumanPreSuspensionData::load_qc(
    selected_omes = "epigen-methylcap-seq",
    epigen = TRUE,
    repo_local_dir = repo_local_dir
  )

  for(tissue in names(data)){
    feature_metadata = data[[tissue]][["epigen-methylcap-seq"]][["qc_norm"]] %>%
      tibble::rownames_to_column("feature_id") %>%
      select(feature_id)

    write_with_path_name(
      actual_data_object = feature_metadata,
      local_path = file.path(repo_local_dir, "data", "tmp", "freeze",
                             "epigenomics", "metadata"),
      tissue = tissue,
      ome = "epigen-methylcap-seq",
      data_category = "metadata",
      data_details = "features",
      version = "1.4",
      return_name_only = FALSE
    )
  }
}
