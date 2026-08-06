#' @title Human Feature to Gene Conversion Table
#'
#' @description A \code{data.table} used to convert from feature IDs to various
#'   identifiers.
#'
#' @usage HUMAN_FEATURE_TO_GENE
#'
#' @format A sorted \code{data.table} with 1,920,618 rows and 10 columns:
#'
#' \describe{
#'   \item{assay}{factor; the assay (ome). One of "epigen-atac-seq", "epigen-methylcap-seq", "metab", "prot-ol", "prot-ph", "prot-pr", or "transcript-rna-seq".}
#'   \item{feature_id}{factor; the feature identifier.}
#'   \item{entrez_gene}{factor; Entrez gene identifier.}
#'   \item{gene_symbol}{factor; gene symbol.}
#'   \item{ensembl_gene}{factor; ensembl gene identifier.}
#'   \item{uniprot}{factor; UniProt identifier.}
#'   \item{refmet_name}{factor; RefMet metabolite name, from the pinned RefMet snapshot.}
#'   \item{refmet_id}{factor; RefMet metabolite identifier, from the pinned RefMet snapshot.}
#'   \item{kegg_id}{factor; KEGG identifier, from the pinned RefMet/KEGG snapshot.}
#'   \item{flanking_sequence}{factor; flanking sequence. Only applicable if assay is "prot-ph".}
#' }
#'
#' @source Built by Stage 1 step 07 of the precovid-repro pipeline. The
#'   \code{refmet_name}, \code{refmet_id} and \code{kegg_id} columns come from a
#'   pinned offline RefMet/KEGG snapshot rather than a live Metabolomics Workbench
#'   query, so the mapping does not move with those databases.
#'
#' @keywords datasets
"HUMAN_FEATURE_TO_GENE"
