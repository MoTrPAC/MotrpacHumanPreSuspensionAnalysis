#' Generate QC-normalized Olink proteomics tables
#'
#' This function generates quality-controlled and batch-corrected normalization
#' tables for Olink targeted proteomics (Prot-OL) data. It downloads raw NPX
#' values and associated metadata from cloud storage, applies feature- and
#' sample-level filtering, removes outliers, and corrects for technical batch
#' effects prior to writing standardized QC-normalized outputs to disk.
#'
#' The resulting QC-normalized matrices are used as inputs for downstream
#' differential analysis and visualization workflows.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory where raw data are
#' downloaded and normalized outputs are written.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Download raw Olink NPX results, sample metadata, and protein metadata
#'   \item Reshape NPX data from wide to long and merge with metadata
#'   \item Filter proteins based on missingness frequency (< 80\%)
#'   \item Pivot data back to wide format and remove samples with missing values
#'   \item Remove predefined sample outliers
#'   \item Restrict analysis to acute study participants
#'   \item Construct covariate matrices using \code{process_covariates}
#'   \item Remove technical batch effects using \code{limma::removeBatchEffect}
#'   \item Write QC-normalized expression matrices and metadata to disk
#' }
#'
#' Batch correction preserves biological design covariates while regressing out
#' technical covariates defined in the covariate specification file.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' QC-normalized Olink proteomics matrices and associated metadata are written
#' to disk.
#'
#' @note
#' Feature filtering thresholds and outlier definitions follow consortium-wide
#' standards and are consistent with prior internal analyses.
#'
#' @keywords internal proteomics normalization
#' @author christopher jin

generate_prot_ol_qc_norm = function(repo_local_dir){
  desired_ome = 'prot-ol'; tissue = 'blood'
  local_path = file.path(repo_local_dir, "data", "tmp/")
  olink_data_path = "gs://motrpac-data-hub/human-precovid/results/proteomics-targeted/t02-plasma/prot-ol/motrpac_human-precovid_t02-plasma_prot-ol_results_v1.0.txt"
  olink_data = MotrpacBicQC::dl_read_gcp(olink_data_path, tmpdir = local_path, check_first = T)
  sample_metadata = MotrpacBicQC::dl_read_gcp('gs://motrpac-data-hub/human-precovid/results/proteomics-targeted/t02-plasma/prot-ol/motrpac_human-precovid_t02-plasma_prot-ol_metadata-samples_v1.0.txt',
                                tmpdir = local_path, check_first = T) %>%
    dplyr::rename(vialLabel = sample_id)
  protein_metadata = MotrpacBicQC::dl_read_gcp('gs://motrpac-data-hub/human-precovid/results/proteomics-targeted/t02-plasma/prot-ol/motrpac_human-precovid_t02-plasma_prot-ol_metadata-proteins_v1.0.txt',
                                 tmpdir = local_path, check_first = T ) %>%
    dplyr::rename(OlinkID = olink_id)

  metadata_path = paste0(local_path, "freeze/proteomics/metadata/")
  qc_norm_path = paste0(local_path, "freeze/proteomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  #this filtering is borrowed from Jacob Barber's code
  message("Generating normalized matrixes for prot-ol")

  t_dat = olink_data %>% t() %>% as.data.frame()
  t_dat = t_dat %>% dplyr::mutate(vialLabel = rownames(t_dat)) %>% dplyr::select(vialLabel, everything()) #reorder vialLabel
  colnames(t_dat) = t_dat[1,] #set colnames then remove the first row
  t_dat = t_dat[-1,]

  t_dat = t_dat %>% rename(vialLabel = olink_id) %>% dplyr::mutate_at(vars(starts_with('OID')), as.numeric)
  raw_olink = t_dat %>% tidyr::pivot_longer(cols=-vialLabel, names_to = 'OlinkID', values_to = 'NPX')
  raw_olink = raw_olink %>% dplyr::left_join(sample_metadata, by='vialLabel')
  raw_olink = raw_olink %>% dplyr::left_join(protein_metadata, by='OlinkID')
  raw_olink = raw_olink %>% dplyr::filter(missing_freq < 0.8)

  raw_olink_wide = raw_olink %>%
    dplyr::select("vialLabel","OlinkID","plate_id","NPX") %>%
    tidyr::pivot_wider(names_from = "OlinkID", values_from = "NPX") %>%
    na.omit()

  #here we remove outliers
  outliers_data = OUTLIERS$vialLabel
  pheno_data_parsed = load_pheno(repo_local_dir, load_acute_only = FALSE)$pheno_data

  raw_olink_wide_manifest <- dplyr::left_join(raw_olink_wide, pheno_data_parsed, by='vialLabel') %>%
    dplyr::select(pid, visitcode, plate_id, vialLabel, BMI, calculatedAge, Timepoint,study, sex_psca,starts_with('OID'))
  norm_adu <- raw_olink_wide_manifest %>%
    dplyr::filter(study=='01') %>%
    dplyr::filter(!vialLabel %in% as.character(outliers_data)) %>%
    dplyr::select(vialLabel, starts_with('OID'))

  new_norm_table <- norm_adu %>% dplyr::select(vialLabel, starts_with('OID')) %>% t()
  colnames(new_norm_table) <- new_norm_table[1,] #set viallabels as colnames
  new_norm_table <- new_norm_table[-1,] %>% as.data.frame() %>%
    dplyr::mutate(across(starts_with('1'), ~as.numeric(.)) )

  sample_metadata_output = sample_metadata[sample_metadata$vialLabel %in% colnames(new_norm_table),]
  protein_metadata_output = protein_metadata %>%
    dplyr::rename(feature_id = OlinkID) %>%
    dplyr::filter(feature_id %in% rownames(new_norm_table))

  wider_metadata = merge(sample_metadata_output, pheno_data_parsed, by = 'vialLabel')
  process_metadata = process_covariates(meta = wider_metadata,
                                        selected_ome = desired_ome,
                                        tissue_input = tissue)
  meta = process_metadata$metadata
  technical_cov = paste(process_metadata[["technical_cov"]]$covariate, collapse = " + ")
  design_cov = paste(process_metadata[["design_cov"]], collapse = " + ")
  message(tissue, ";", desired_ome, ";technical: ", technical_cov, ";design: ", design_cov)

  new_norm_table = limma::removeBatchEffect(new_norm_table,
                                     covariates = model.matrix(as.formula(paste("~ ", technical_cov)), data = meta),
                                     design = model.matrix(as.formula(paste("~ ", design_cov)), data = meta))

  oids <- rownames(new_norm_table)
  new_norm_table <- as.data.frame(new_norm_table)
  new_norm_table$feature_id <- oids
  new_norm_table <- new_norm_table %>% dplyr::select(feature_id, everything()) #set feature_id as first col

  write_with_path_name(sample_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
  write_with_path_name(new_norm_table, local_path = qc_norm_path, ome = desired_ome, tissue = tissue, data_category = 'qc-norm', data_details = 'log2')
  write_with_path_name(protein_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features')
}

