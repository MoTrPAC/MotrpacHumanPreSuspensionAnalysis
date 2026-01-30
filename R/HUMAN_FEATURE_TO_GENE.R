#' @title Human Feature to Gene Conversion Table
#'
#' @description A \code{data.table} used to convert from feature IDs to various
#'   identifiers.
#'
#' @usage HUMAN_FEATURE_TO_GENE
#'
#' @format A sorted \code{data.table} with 1,958,568 rows and 11 columns:
#'
#' \describe{
#'   \item{assay}{factor; the assay (ome). One of "epigen-atac-seq",
#'   "epigen-methylcap-seq", "metab", "prot-ol", "prot-ph", "prot-pr", or
#'   "transcript-rna-seq".}
#'   \item{feature_id}{factor; the feature identifier.}
#'   \item{entrez_gene}{factor; Entrez gene identifier.}
#'   \item{gene_symbol}{factor; gene symbol.}
#'   \item{ensembl_gene}{factor; ensembl gene identifier.}
#'   \item{custom_annotation}{factor; custom feature annotation. Only applicable
#'   if assay is "epigen-atac-seq" or "epigen-methylcap-seq".}
#'   \item{relationship_to_gene}{numeric; only applicable if assay is
#'   "epigen-atac-seq" or "epigen-methylcap-seq".}
#'   \item{uniprot}{factor; UniProt identifier.}
#'   \item{refmet_name}{factor; RefMet metabolite identifier.}
#'   \item{kegg_id}{factor; KEGG identifier.}
#'   \item{flanking_sequence}{factor; flanking sequence. Only applicable if
#'   assay is "prot-ph".}
#' }
#'
#' @keywords datasets
"HUMAN_FEATURE_TO_GENE"
