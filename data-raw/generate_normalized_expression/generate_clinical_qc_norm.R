
#' Generate QC-normalized clinical chemistry tables
#'
#' This function differs from the other \code{generate_*_qc_norm()} functions in
#' this package. For molecular omes (transcriptomics, proteomics, metabolomics,
#' etc.), raw data are pulled from cloud storage and processed here. Clinical
#' chemistry data are instead curated and normalized upstream by the MoTrPAC
#' clinical group; this function ingests those already-processed objects directly
#' from \code{MotrpacHumanPreSuspensionData} rather than recomputing normalization
#' from raw values.
#'
#' The only real functionality here is to make a feature-metadata file to accompany this raw data.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory where outputs are
#' written.
#'
#' @return
#' This function is called for its side effects and returns no object.
#'
#' @keywords internal clinical normalization
#' @author christopher jin

generate_clinical_qc_norm = function(config_file = "~/config.json"){

  if(!file.exists(config_file))
    stop(config_file, " not found. This function requires a config file to locate precovid_repo_path.")

  config = jsonlite::fromJSON(config_file)
  repo_local_dir = config$precovid_repo_path
  metadata_output = file.path(repo_local_dir, "data", "tmp", "freeze", "clinical_chemistry", "metadata")
  qc_output = file.path(repo_local_dir, "data", "tmp", "freeze", "clinical_chemistry", "qc-norm")

  dir.create(metadata_output, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_output, recursive = TRUE, showWarnings = FALSE)

  require(MotrpacHumanPreSuspensionData)
  require(dplyr)
  require(here)

  clin_chemistry = MotrpacHumanPreSuspensionData::load_clinical_data()[["chemistry"]]

  #note these are all in absolute terms, we log2 transform them and then fit linear models on them pretty directly.

  #there is one missing value though. we treat it as missing completely at random.
  log_clin_chem_df = clin_chemistry %>%
    dplyr::bind_rows() %>%
    tibble::column_to_rownames("analyte_name") %>%
    dplyr::mutate(across(everything(), ~log2(.x)))

  # log-QC-norm
  log_cln_out = log_clin_chem_df %>% tibble::rownames_to_column("feature_id") %>%
    select(feature_id, everything())

  # note: .write_with_path_names doesn't really work here because the clinical chemistry is separate and there
  # isn't really a tissue_code corresponding that we can automatically apply.
  file_type = ".txt"
  all_file_header = "human-precovid-sed-adu" #this is the base structure for all files within the phase.
  tissue_code = "t02-plasma"

  file_name = paste(all_file_header, tissue_code, "clinical_chemistry", "qc-norm", "log2-transformed", sep = "_")
  file_name = paste0(qc_output, "/", file_name, "_v", "1.3", file_type)
  write.table(log_cln_out, file = file_name, row.names = F, sep = '\t', quote = F)

  # as we get ready to annotate this, there are some clinical features that are proteins
  prot_clinical_feautres = c("Insulin", "Glucagon", "CK")
  metab_clinical_features = setdiff(log_cln_out$feature_id, prot_clinical_feautres)

  prot_annotation = .annotate_clinical_prot(prot_clinical_feautres) %>%
    dplyr::mutate(assay = "prot-clinical")
  metab_annotation = .annotate_clinical_metab(metab_clinical_features) %>%
    dplyr::mutate(assay = "metab-t-clinical")
  all_annotation = dplyr::bind_rows(prot_annotation, metab_annotation)

  feature_metadata = log_cln_out %>%
    dplyr::select(feature_id) %>%
    dplyr::left_join(all_annotation, by = "feature_id")

  feature_metadata_file = paste(all_file_header, tissue_code, "clinical_chemistry", "metadata", "features", sep = "_")
  feature_metadata_file = paste0(metadata_output, "/", feature_metadata_file, "_v", "1.3", file_type)
  write.table(feature_metadata, file = feature_metadata_file, row.names = F, sep = '\t', quote = F)
}

.annotate_clinical_prot = function(prot_analyte_names){
  #analyte names don't match Ensembl external_gene_name; map manually to gene symbols first
  prot_analyte_to_gene = c("Insulin" = "INS", "Glucagon" = "GCG", "CK" = "CKM")

  ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  biomart_anno = biomaRt::getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprotswissprot"),
    filters = "external_gene_name",
    values = prot_analyte_to_gene[prot_analyte_names],
    mart = ensembl
  ) %>%
    dplyr::distinct() %>%
    dplyr::filter(uniprotswissprot != "") %>%
    dplyr::mutate(ensembl_gene = str_remove(ensembl_gene_id, "\\..*"),
                  entrez_gene = as.character(entrezgene_id),
                  gene_symbol = external_gene_name,
                  feature_id = names(prot_analyte_to_gene)[match(external_gene_name, prot_analyte_to_gene)]) %>%
    dplyr::select(feature_id, gene_symbol, ensembl_gene, entrez_gene, uniprot = uniprotswissprot)

  return(biomart_anno)
}

.annotate_clinical_metab = function(metab_analyte_names){
  mets = paste(metab_analyte_names, collapse = "\n")

  h = curl::new_handle()
  curl::handle_setform(h, metabolite_name = mets)
  req = curl::curl_fetch_memory(
    "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_minID.php",
    handle = h
  )

  refmet_result = utils::read.table(
    text = rawToChar(req$content),
    header = TRUE,
    na.strings = "-",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    sep = "\t"
  )

  refmet_result[is.na(refmet_result)] = '-'
  refmet_result[refmet_result == ''] = '-'

  metab_df = data.frame(feature_id = metab_analyte_names,
                        lookup_refmet = metab_analyte_names)

  refmet_output = refmet_result %>%
    dplyr::transmute(
      lookup_refmet = Input.name,
      refmet_name_std = Standardized.name,
      refmet_id = RefMet_ID,
      kegg_id = KEGG_ID
    ) %>%
    dplyr::left_join(metab_df, by = "lookup_refmet", relationship = "many-to-many") %>%
    dplyr::transmute(
      feature_id,
      refmet_name = dplyr::na_if(refmet_name_std, "-"),
      refmet_id = dplyr::na_if(refmet_id, "-"),
      kegg_id = dplyr::na_if(kegg_id, "-")
    ) %>%
    dplyr::distinct()

  return(refmet_output)
}
