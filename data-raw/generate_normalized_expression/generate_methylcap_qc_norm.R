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

generate_methylcap_qc_norm = function(config_file = "~/config.json"){
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
      dplyr::select(feature_id)
    methyl_annotation = .annotate_methyl_features(feature_metadata)

    write_with_path_name(
      actual_data_object = methyl_annotation,
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

# get_peak_annotations_hs is run per chromosome to stay within memory limits —
# passing the full feature set at once exhausts RAM on large methylcap datasets.
.annotate_methyl_features = function(feature_metadata){
  feature_input = feature_metadata %>%
    dplyr::mutate(chrom = gsub(":.*", "", feature_id),
                  start = as.numeric(gsub(".*:|-.*", "", feature_id)),
                  end = as.numeric(gsub(".*-", "", feature_id))) %>%
    dplyr::mutate(chrom = gsub("chrM", "chrMT", chrom)) %>%
    data.table::as.data.table()

  methyl_peakdf = pre_cawg_get_peak_annotations_hs(feature_input)
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 105)

  attributes <- c("ensembl_gene_id", "entrezgene_id", "external_gene_name")
  methyl_lookup_df <- biomaRt::getBM(attributes = attributes,
                                     filters = "ensembl_gene_id",
                                     values = methyl_peakdf$ensembl_gene,
                                     mart = ensembl) %>%
    dplyr::full_join(methyl_peakdf, c("ensembl_gene_id" = "ensembl_gene")) %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(external_gene_name == "", NA, external_gene_name)) %>%
    dplyr::select(feature_id,
                  gene_symbol,
                  ensembl_gene = ensembl_gene_id,
                  entrez_gene = entrezgene_id,
                  custom_annotation,
                  relationship_to_gene) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::slice(match(min(entrez_gene), entrez_gene)) %>%
    dplyr::mutate(entrez_gene = as.character(entrez_gene),
                  platform = "epigen-methylcap-seq")

  return(methyl_lookup_df)
}




