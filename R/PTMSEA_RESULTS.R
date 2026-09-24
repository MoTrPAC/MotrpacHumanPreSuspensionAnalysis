#' @title PTM-SEA Results
#'
#' @description PTM-SEA enrichment of PTMsigDB signatures in the prot-ph
#'   differential-analysis z-statistics (see `PTMSEA_INPUT`), run with the Broad
#'   Institute PTM-SEA / ssGSEA2.0 container, for the acute-exercise EE-CON and
#'   RE-CON contrasts. Muscle carries three post-exercise timepoints crossed
#'   with the two exercise modalities; adipose was sampled at `post_3.5_4_hr`
#'   only. Same layout as \code{\link{CAMERA_RESULTS}}, with `NES` in place of
#'   `t`, `df` and `z.std`, so it can be passed to
#'   \code{\link{plot_enrich_heatmap}}.
#'
#' @usage
#' PTMSEA_RESULTS
#'
#' @format A \code{data.frame} with 3,910 rows (muscle: 506 signatures x 6
#'   contrasts; adipose: 437 signatures x 2 contrasts) and 17 columns:
#'   \describe{
#'     \item{\code{tissue}}{factor; `"muscle"` or `"adipose"`.}
#'     \item{\code{assay}}{factor; `"prot-ph"`.}
#'     \item{\code{contrast_type}}{factor; `"exercise_with_controls"`.}
#'     \item{\code{contrast}}{factor; the contrast, as in
#'       \code{CONTRAST_CONVERTER$contrast}.}
#'     \item{\code{contrast_short}}{factor; abbreviated contrast label.}
#'     \item{\code{collection}}{factor; `"PTMSIGDB"`.}
#'     \item{\code{database}}{factor; the PTMsigDB signature category, the
#'       prefix of \code{set} (e.g. `"KINASE-PSP"`, `"PERT-P100-DIA2"`).}
#'     \item{\code{set_id}}{factor; the PTMsigDB signature ID. Not an ID from
#'       \code{\link{SET_TO_ID}}.}
#'     \item{\code{set}}{factor; the PTMsigDB signature ID.}
#'     \item{\code{set_short}}{factor; signature name followed by
#'       \code{database} in parentheses.}
#'     \item{\code{set_size}}{integer; number of signature sites measured.}
#'     \item{\code{set_size_DB}}{integer; number of sites in the signature, as
#'       reported by PTM-SEA.}
#'     \item{\code{size_ratio}}{numeric; \code{set_size / set_size_DB}.}
#'     \item{\code{direction}}{factor; `"Up"` if \code{NES > 0}, else
#'       `"Down"`.}
#'     \item{\code{NES}}{numeric; normalized enrichment score.}
#'     \item{\code{p_value}}{numeric; nominal p-value.}
#'     \item{\code{adj_p_value}}{numeric; PTM-SEA FDR, adjusted within tissue
#'       and contrast.}
#'   }
#'
#' @source Muscle: precovid-analyses `figures/muscle/Figure5/` (Natalie Clark).
#'   Adipose: Cheehoon Ahn. Provenance in `data-raw/PTMSEA/README.md`.
#'
#' @seealso \code{\link{CAMERA_RESULTS}}, \code{\link{plot_enrich_heatmap}},
#'   \code{\link{PTMSEA_INPUT}}
#'
#' @keywords datasets proteomics
"PTMSEA_RESULTS"
