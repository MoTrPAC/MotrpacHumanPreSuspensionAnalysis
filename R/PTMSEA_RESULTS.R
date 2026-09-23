#' @title PTM-SEA Results
#'
#' @description PTM-SEA enrichment of PTMsigDB signatures in the prot-ph
#'   differential-analysis z-statistics (see `PTMSEA_INPUT`), run with the Broad
#'   Institute PTM-SEA / ssGSEA2.0 container. One row per signature, one column
#'   per acute-exercise contrast (EE-CON and RE-CON). Muscle carries three
#'   post-exercise timepoints crossed with the two exercise modalities; adipose
#'   was sampled at `post_3.5_4_hr` only.
#'
#'   Each element is the combined PTM-SEA output parsed into a `GCT` object, the
#'   S4 class defined by the `cmapR` package. This package does not depend on
#'   `cmapR`: the object loads without it, but the `cmapR` accessors need
#'   `cmapR` attached. The slots are reachable as `@mat`, `@rdesc`, `@cdesc`,
#'   `@rid` and `@cid` either way.
#'
#' @usage
#' PTMSEA_RESULTS
#'
#' @format A named list of two `cmapR` `GCT` objects, `"muscle"` (506
#'   signatures x 6 contrasts) and `"adipose"` (437 signatures x 2 contrasts).
#'   Row IDs are PTMsigDB signature IDs. The `mat` slot holds the normalized
#'   enrichment scores (NES). The `rdesc` slot holds, per signature:
#'   \describe{
#'     \item{Signature.set.description}{Signature sites and supporting PubMed IDs}
#'     \item{Signature.set.size}{Number of sites in the signature}
#'     \item{Signature.set.overlap.percent.<contrast>}{Percent of signature sites measured}
#'     \item{Signature.set.overlap.<contrast>}{Measured signature sites, `|`-separated}
#'     \item{No.columns.scored}{Number of contrasts scored}
#'     \item{pvalue.<contrast>}{Nominal p-value}
#'     \item{fdr.pvalue.<contrast>}{FDR-adjusted p-value}
#'   }
#'
#' @source Muscle: precovid-analyses `figures/muscle/Figure5/` (Natalie Clark).
#'   Adipose: Cheehoon Ahn. Provenance in `data-raw/PTMSEA/README.md`.
#'
#' @keywords datasets proteomics
"PTMSEA_RESULTS"
