#' @title Human Feature to Gene Conversion Table
#'
#' @description A \code{data.table} used to convert from feature IDs to various
#'   identifiers.
#'
#' @usage HUMAN_FEATURE_TO_GENE
#'
#' @format A sorted \code{data.table} with 1,920,618 rows and 12 columns:
#'
#' \describe{
#'   \item{assay}{character; the assay (ome). One of "epigen-atac-seq", "epigen-methylcap-seq", "metab", "prot-clinical", "prot-ol", "prot-ph", "prot-pr", or "transcript-rna-seq".}
#'   \item{feature_id}{factor; the feature identifier.}
#'   \item{entrez_gene}{factor; Entrez gene identifier.}
#'   \item{gene_symbol}{factor; gene symbol.}
#'   \item{ensembl_gene}{factor; ensembl gene identifier.}
#'   \item{uniprot}{factor; UniProt identifier.}
#'   \item{refmet_name}{factor; RefMet metabolite name, from the pinned RefMet snapshot.}
#'   \item{refmet_id}{factor; RefMet metabolite identifier, from the pinned RefMet snapshot.}
#'   \item{kegg_id}{factor; KEGG identifier, from the pinned RefMet/KEGG snapshot.}
#'   \item{custom_annotation}{factor; the region of the assigned gene the peak falls in — one of "Promoter (<=1kb)", "Promoter (1-2kb)", "5' UTR", "3' UTR", "Exon", "Intron", "Overlaps Gene", "Upstream (<5kb)", "Downstream (<5kb)" or "Distal Intergenic". \code{NA} unless assay is "epigen-atac-seq" or "epigen-methylcap-seq".}
#'   \item{relationship_to_gene}{numeric; the signed distance in base pairs from the peak to the assigned gene, \code{0} where the peak overlaps it. \code{NA} unless assay is "epigen-atac-seq" or "epigen-methylcap-seq".}
#'   \item{flanking_sequence}{factor; flanking sequence. Only applicable if assay is "prot-ph".}
#' }
#'
#' @details \code{custom_annotation} and \code{relationship_to_gene} describe where an
#'   ATAC-seq or MethylCap-seq peak sits relative to the gene it was assigned to. Without
#'   them a peak in a promoter and a peak 40 kb into an intron are indistinguishable once
#'   mapped, since both carry only the gene. Both are derived from the peak coordinates in
#'   the \code{feature_id} rather than measured per tissue, so they take the same value in
#'   every tissue a peak appears in.
#'
#'   Per-tissue measurements are deliberately not carried here. This table is keyed on
#'   \code{(assay, feature_id)} with no tissue column, so it could only hold a collapse
#'   across tissues. The phosphosite localization flag (\code{confident_site}) is the
#'   case in point: it was added in 2.0.3 and removed in 2.0.7 because muscle and adipose
#'   disagree on 859 of their 7,865 shared prot-ph sites. Read it per tissue from
#'   \code{*_PROT_PH_QC$feature_metadata} in \pkg{MotrpacHumanPreSuspensionData}.
#'
#' @source Built by Stage 1 step 07 of the motrpac-human-presuspension-repro pipeline. The
#'   \code{refmet_name}, \code{refmet_id} and \code{kegg_id} columns come from a
#'   pinned offline RefMet/KEGG snapshot rather than a live Metabolomics Workbench
#'   query, so the mapping does not move with those databases. The peak annotations in
#'   \code{custom_annotation} and \code{relationship_to_gene} are produced by
#'   \code{ChIPseeker} against the pinned Ensembl v105 \code{TxDb}.
#'
#' @keywords datasets
"HUMAN_FEATURE_TO_GENE"
