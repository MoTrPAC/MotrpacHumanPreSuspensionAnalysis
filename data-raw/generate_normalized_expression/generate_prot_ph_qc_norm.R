#' Generate QC-normalized phosphoproteomics (Prot-PH) tables
#'
#' This function generates quality-controlled and normalized phosphoproteomics
#' (Prot-PH) data matrices for downstream analysis. The workflow was originally
#' derived from the untargeted phosphoproteomics pipeline and includes multiple
#' stages of filtering, normalization, replicate handling, and batch correction.
#'
#' Raw phosphoproteomics ratio results are loaded from cloud storage, filtered to
#' remove contaminants and duplicated PTM identifiers, normalized using a
#' median-based strategy, and corrected for technical batch effects. Tissue-
#' specific post-processing steps are applied to appropriately handle replicate
#' measurements.
#'
#' @param repo_local_dir
#' Character scalar specifying the local repository directory where raw data are
#' downloaded and QC-normalized outputs are written.
#'
#' @details
#' For each supported tissue, the function performs the following steps:
#' \enumerate{
#'   \item Download raw phosphoproteomics ratio results and vial-level metadata
#'   \item Remove duplicated PTM identifiers and known contaminants
#'   \item Construct a GCT object linking expression data with row and column metadata
#'   \item Restrict analysis to acute study participants
#'   \item Apply median normalization to phosphoproteomics ratios
#'   \item Remove predefined sample outliers
#'   \item Define technical and biological covariates via \code{process_covariates}
#'   \item Remove technical batch effects while preserving biological design factors
#'   \item Collapse intra-site and inter-site technical replicates (muscle only)
#'   \item Filter features based on missingness thresholds
#'   \item Write QC-normalized matrices and metadata to disk
#' }
#'
#' For muscle tissue, replicate measurements are averaged according to predefined
#' intra-site and inter-site conventions prior to missingness filtering.
#'
#' @return
#' This function is called for its side effects and returns no object.
#' QC-normalized phosphoproteomics matrices and associated metadata are written
#' to disk.
#'
#' @note
#' This function requires \code{cmapR} at runtime for handling GCT objects.
#' The package is checked dynamically and is not declared as a formal dependency.
#'
#' @keywords internal phosphoproteomics normalization
#' @author christopher jin

generate_prot_ph_qc_norm = function(repo_local_dir){
  local_path = repo_local_dir

  desired_ome = 'prot-ph'; tissue_types = c('muscle', 'adipose')
  ome_vial_meta = list()
  ome_vial_meta[["muscle"]] = "gs://motrpac-data-hub/human-precovid/results/proteomics-untargeted/t10-muscle/prot-ph/motrpac_human-precovid_t10-muscle_prot-ph_vial-metadata_v1.0-pnbi.txt"
  ome_vial_meta[["adipose"]] = "gs://motrpac-data-hub/human-precovid/results/proteomics-untargeted/t07-adipose/prot-ph/motrpac_human-precovid_t07-adipose_prot-ph_vial-metadata_v1.0.txt"

  pheno_data_parsed = load_pheno(load_acute_only = FALSE)$pheno_data
  metadata_path = paste0(local_path, "freeze/proteomics/metadata/")
  qc_norm_path = paste0(local_path, "freeze/proteomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  for (tissue in tissue_types){
    message(paste("Generating normalized matrixes for phospho-ph", tissue))
    file_load = .find_path_name(desired_ome, tissue = tissue, data_type = 'ratio-results', version = "1.0") #find version 1 of the ratio results
    phospho = MotrpacBicQC::dl_read_gcp(file_load, tmpdir = local_path)
    phospho = phospho[!duplicated(phospho$ptm_id),] #remove duplicates
    phospho = phospho %>% dplyr::filter(phospho$is_contaminant == "FALSE")

    mat = phospho %>% dplyr::select(matches("^[0-9]")) %>% as.data.frame()
    rdesc <- phospho %>% dplyr::select(!matches("^[0-9]"))
    rdesc <- rdesc %>% dplyr::select(!starts_with("pool"))

    rdesc$id <- rdesc$ptm_id
    rownames(mat) = rdesc$ptm_id

    tmt_metadata = MotrpacBicQC::dl_read_gcp(ome_vial_meta[[tissue]], tmpdir = local_path) %>%
      dplyr::rename(vialLabel= vial_label) %>%
      dplyr::mutate(id = vialLabel) %>%
      dplyr::mutate(vialLabel = gsub("\\.1", "", vialLabel))
    meta_merge = dplyr::left_join(tmt_metadata, pheno_data_parsed, by = 'vialLabel')
    meta_merge <- meta_merge[match(colnames(mat),meta_merge$id),]

    if(ncol(mat)!=nrow(meta_merge) | nrow(mat)!=nrow(rdesc)){
      stop("Column or row annotations do not match
		     the matrix size.")}

    #this is non-batch effect corrected, non-normalized dataset
    phospho_nonnorm <- cmapR::GCT(mat=as.matrix(mat),
                                  rdesc=rdesc,
                                  cdesc=meta_merge,
                                  rid =rownames(mat),
                                  cid = colnames(mat))
    phospho_nonnorm <- cmapR::subset_gct(phospho_nonnorm, cid = which(phospho_nonnorm@cdesc$study == "01")) #remove HA participants
    #This is non-batch effect corrected, median-normalized dataset
    phospho_mednorm <- .median_mad_norm(phospho_nonnorm, mad = F)

    remove = OUTLIERS$vialLabel
    phospho_mednorm_no_outliers <- cmapR::subset_gct(phospho_mednorm, cid=which(!(phospho_mednorm@cdesc$vialLabel %in% remove)))

    #batch correction
    phospho_sample_meta = phospho_mednorm_no_outliers@cdesc
    process_metadata = process_covariates(meta = phospho_sample_meta,
                                          selected_ome = desired_ome,
                                          tissue_input = tissue)
    meta = process_metadata$metadata
    technical_cov = paste(process_metadata[["technical_cov"]]$covariate, collapse = " + ")
    design_cov = paste(process_metadata[["design_cov"]], collapse = " + ")
    message(tissue, ";", desired_ome, ";technical: ", technical_cov, ";design: ", design_cov)

    phospho_mednorm_no_outliers@mat = limma::removeBatchEffect(phospho_mednorm_no_outliers@mat,
                                                               covariates = model.matrix(as.formula(paste("~ ", technical_cov)), data = meta),
                                                               design = model.matrix(as.formula(paste("~ ", design_cov)), data = meta))
    if(tissue == 'muscle'){
      matr <- phospho_mednorm_no_outliers@mat %>% as.data.frame()
      #Per patrick: INTRA site sites end in .1 (and the other version of the rep will have the same first few digits)
      #INTER site replicates end in 7 (and the other version of the rep will have the same first few digits)
      intrasite_replicates = phospho_mednorm_no_outliers@cdesc$id[grepl(".1$", phospho_mednorm_no_outliers@cdesc$id)]
      for(i in intrasite_replicates){
        j <- gsub("\\.1$", "", i)
        matr[i] <- rowMeans(matr[, c(i, j)], na.rm = TRUE)
      }
      #averaging replicates across chemical sites
      intersite_replicates = phospho_mednorm_no_outliers@cdesc$id[grepl("7$", phospho_mednorm_no_outliers@cdesc$id)]
      for(j in intersite_replicates){
        i <- paste0(stringr::str_sub(j, end = -2), "2", sep = "")
        matr[i] <- rowMeans(matr[, c(i, j)], na.rm = TRUE)
      }
      #removing columns which have been averaged
      matr <- dplyr::select(matr, -contains(".1"))
      matr <- dplyr::select(matr, -all_of(intersite_replicates))
      matr <- as.matrix(matr)

      cdescc = phospho_sample_meta %>%
        dplyr::filter(!grepl("\\.1$", id)) %>%
        dplyr::filter(id %in% colnames(matr))
      cdescc <- cdescc[match(colnames(matr),cdescc$id),]

      if(ncol(matr)!=nrow(cdescc)|nrow(mat)!=nrow(rdesc)){
        stop("Column or row annotations do not match
		     the matrix size. (after avg)")}

      phospho_ave <- cmapR::GCT(mat=matr, rdesc= data.frame(rdesc), cdesc=data.frame(cdescc), rid =rownames(matr), cid = colnames(matr))
      PH <- phospho_ave %>% .remove_na(0.7)
    }else{
      PH = phospho_mednorm_no_outliers %>% .remove_na(0.7)
    }
    phospho_output <- as.data.frame(PH@mat)
    phospho_output$feature_id = rownames(phospho_output)
    phospho_output = phospho_output %>% dplyr::select(feature_id, everything())

    sample_metadata_output = PH@cdesc %>% dplyr::select(c("vialLabel", "tmt_plex", "tmt16_channel"))

    feature_metadata_output = PH@rdesc %>%
      dplyr::rename(feature_id = ptm_id) %>%
      dplyr::select(feature_id, everything()) %>%
      dplyr::filter(feature_id %in% rownames(phospho_output)) %>%
      .annotate_prot_ph()

    write_with_path_name(sample_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
    write_with_path_name(phospho_output, local_path = qc_norm_path, ome = desired_ome, tissue = tissue, data_category = 'qc-norm', data_details = 'log2-mn')
    write_with_path_name(feature_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features', version = "1.4")
  }
}

.annotate_prot_ph = function(feature_metadata_output){

  #feature_id is the ptm_id; protein_id is the UniProt accession used for lookup
  prot_ph = feature_metadata_output %>%
    dplyr::mutate(platform = "prot-ph",
                  uniprot = protein_id)

  # human UniProt database download
  uniprot_db = .get_uniprot_mapping()

  ####compile prot-ph and annotate####
  prot_features = prot_ph %>%
    dplyr::mutate(uniprot_lookup = str_remove(uniprot, "-.*")) %>%
    dplyr::left_join(uniprot_db, by = c("uniprot_lookup" = "UniProtKB-AC")) %>%
    dplyr::mutate(ensg_lookup = str_remove(Ensembl, "\\..*"))

  ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprotswissprot")

  prot_lookup_df = biomaRt::getBM(attributes = attributes,
                                   filters = "uniprotswissprot",
                                   values = prot_features$uniprot_lookup,
                                   mart = ensembl) %>%
    dplyr::distinct() %>%
    dplyr::full_join(prot_features, by = c("uniprotswissprot" = "uniprot_lookup")) %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(external_gene_name == "", NA, external_gene_name)) %>%
    dplyr::mutate(ensembl_gene = dplyr::if_else(is.na(ensembl_gene_id),
                                                 str_remove(Ensembl, "\\..*"),
                                                 ensembl_gene_id)) %>%
    dplyr::mutate(entrez_gene = dplyr::if_else(is.na(entrezgene_id),
                                               `GeneID (EntrezGene)`,
                                               as.character(entrezgene_id))) %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol),
                                               str_remove(`UniProtKB-ID`, "_HUMAN"),
                                               gene_symbol)) %>%
    dplyr::select(entrez_gene, feature_id, gene_symbol, uniprot, ensembl_gene, platform) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::slice(match(min(entrez_gene), entrez_gene))

  #these are a few incorrect gene names that are manually corrected
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "SHAN3"] <- "SHANK3"
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "HECD4"] <- "HECTD4"
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "WASH6"] <- "WASH6P"

  prot_lookup_df = prot_lookup_df %>%
    dplyr::select(assay = platform, feature_id, entrez_gene, gene_symbol, ensembl_gene, uniprot)

  return(prot_lookup_df)
}
