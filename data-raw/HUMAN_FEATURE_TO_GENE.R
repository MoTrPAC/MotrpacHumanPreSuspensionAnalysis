library(dplyr)
library(data.table)
library(MotrpacBicQC)
library(here)

#updated following the change to version 1.4 in the feature to gene map.
#see the v1.4 patch notes report for more details

#just for 'write_with_path_name'
#location for output of file
config_file = jsonlite::fromJSON("~/config.json")
PRECOVID_REPO_PATH = config_file$precovid_repo_path

reference_bucket = "gs://pre-cawg/staging_20260511"
gsutil_command = "gsutil" #only if you have gsutil in your $PATH$. Otherwise whereever you put it
all_files = system(command = paste(gsutil_command, "ls -R", reference_bucket), intern = TRUE)
feature_metadata_files = all_files[grep("metadata_features", all_files)]

all_feat_to_gene = list()

for(file in feature_metadata_files){
  current_file = MotrpacBicQC::dl_read_gcp(path = file,
                                           tmpdir = tempdir()) %>%
    mutate(across(everything(), as.character))

  all_feat_to_gene[[file]] = current_file
}

HUMAN_FEATURE_TO_GENE = dplyr::bind_rows(all_feat_to_gene) %>%
  dplyr::arrange(assay, feature_id) %>%
  dplyr::filter(is_named != "FALSE"| is.na(is_named)) %>%
  dplyr::select(assay, feature_id, entrez_gene, gene_symbol, ensembl_gene, uniprot, refmet_name, refmet_id, kegg_id) %>%
  distinct()

setDT(HUMAN_FEATURE_TO_GENE)

cols <- colnames(HUMAN_FEATURE_TO_GENE)
char_cols <- vapply(cols, class, character(1L)) == "character"
char_cols <- setdiff(names(char_cols),
                     c("relationship_to_gene", "assay"))
HUMAN_FEATURE_TO_GENE[, (char_cols) := lapply(.SD, as.factor),
                      .SDcols = char_cols]

setcolorder(x = HUMAN_FEATURE_TO_GENE,
            neworder = c("assay", "feature_id"))

setkeyv(x = HUMAN_FEATURE_TO_GENE,
        cols = c("assay", "feature_id"))

# write.table(HUMAN_FEATURE_TO_GENE,
#             file = file.path(PRECOVID_REPO_PATH, "data", "tmp", "freeze", "resources", "motrpac-mappings-human-feature-to-gene_v1.4.txt"),
#             sep = "\t",
#             row.names = FALSE)

# Save
usethis::use_data(HUMAN_FEATURE_TO_GENE, overwrite = TRUE,
                  compress = TRUE, version = 3)
