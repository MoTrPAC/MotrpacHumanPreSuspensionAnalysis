# The major data objects for this package that didn't neatly fit into one of the other files.

#' Color scheme for sex-based differentiation
#'
#' Named color vector used to consistently encode biological sex in figures.
#' Male and female labels (case-insensitive) are mapped to fixed colors to
#' ensure visual consistency across plots.
#'
#' @usage HUMAN_SEX_COLORS
#' @format A named character vector of hexadecimal color codes.
"HUMAN_SEX_COLORS"


#' Color scheme for tissue-based differentiation
#'
#' Named color vector assigning consistent colors to tissues commonly analyzed
#' in human exercise studies (e.g., blood, skeletal muscle, adipose).
#'
#' @usage HUMAN_TISSUE_COLORS
#' @format A named character vector of hexadecimal color codes.
"HUMAN_TISSUE_COLORS"


#' Tissue abbreviations for compact figure annotation
#'
#' Named character vector providing short tissue abbreviations for use in
#' heatmap annotations, labels, and compact visual encodings.
#'
#' @usage HUMAN_TISSUE_ABBR
#' @format A named character vector of tissue abbreviations.
"HUMAN_TISSUE_ABBR"


#' Color scheme for exercise group differentiation
#'
#' Named color vector used to distinguish exercise intervention groups,
#' including resistance training, endurance training, and controls.
#'
#' @usage HUMAN_EXERCISE_GROUP_COLORS
#' @format A named character vector of hexadecimal color codes.
"HUMAN_EXERCISE_GROUP_COLORS"


#' Color scheme for acute exercise timepoints
#'
#' Named color vector encoding temporal phases of acute exercise experiments.
#' Multiple synonymous labels are intentionally mapped to the same color to
#' emphasize biological phase rather than exact sampling time.
#'
#' @usage HUMAN_ACUTE_TIMEPOINT_COLORS
#' @format A named character vector of hexadecimal color codes.
"HUMAN_ACUTE_TIMEPOINT_COLORS"


#' Color scheme for molecular modality (omics layer) differentiation
#'
#' Named color vector assigning consistent colors to molecular modalities
#' (e.g., transcriptomics, proteomics, metabolomics, epigenomics). Assay-
#' specific identifiers corresponding to the same modality share a color to
#' facilitate cross-platform comparison.
#'
#' @usage HUMAN_OME_COLORS
#' @format A named character vector of hexadecimal color codes.
"HUMAN_OME_COLORS"



#' @title Fuzzy c-means clustering and enrichment analysis results
#'
#' @description
#' Results derived from fuzzy c-means (FCM) clustering of molecular features,
#' together with downstream gene set enrichment analyses. These datasets
#' summarize clustering assignments and cluster-level functional enrichment
#' using both CAMERA-based competitive gene set testing and classical
#' over-representation analysis (ORA).
#'
#' The FCM results provide soft cluster membership estimates, allowing features
#' to partially belong to multiple clusters, while the enrichment results
#' characterize the biological functions associated with each cluster.
#'
#' @usage
#' FCM_CLUSTERS
#' FCM_CAMERA
#' FCM_ORA
#'
#' @format
#' Objects of class \code{data.frame}.
#'
#' \describe{
#'   \item{\code{FCM_CLUSTERS}}{
#'     Feature-level fuzzy c-means clustering results. Typically includes
#'     feature identifiers, cluster centers, and membership values indicating
#'     the degree of association between each feature and each cluster.
#'   }
#'   \item{\code{FCM_CAMERA}}{
#'     Cluster-level functional enrichment results obtained using CAMERA.
#'     Contains statistics such as enrichment direction, test statistics,
#'     p-values, and multiple-testing–adjusted p-values for each gene set
#'     evaluated within each cluster.
#'   }
#'   \item{\code{FCM_ORA}}{
#'     Cluster-level over-representation analysis results. Summarizes gene sets
#'     that are significantly over-represented among features assigned to each
#'     cluster, along with corresponding enrichment statistics.
#'   }
#' }
#'
#' @details
#' These datasets are intended to be used jointly. Fuzzy c-means clustering
#' defines coherent temporal or condition-specific feature patterns, while
#' CAMERA and ORA provide complementary perspectives on the functional
#' interpretation of those clusters.
#'
#' @keywords datasets clustering enrichment
#'
#' @name FCM_RESULTS

#' @rdname FCM_RESULTS
#' @format NULL
#' @usage NULL
"FCM_CLUSTERS"

#' @rdname FCM_RESULTS
#' @format NULL
#' @usage NULL
"FCM_CAMERA"

#' @rdname FCM_RESULTS
#' @format NULL
#' @usage NULL
"FCM_ORA"



#' @title Molecular pathway signatures
#'
#' @description
#' Collection of molecular pathway and gene set signatures derived from
#' MSigDB. This dataset aggregates curated pathway definitions spanning
#' multiple biological databases and is used as the reference universe for
#' gene set–based analyses throughout the package.
#'
#' Molecular signatures are provided in a standardized format to support
#' enrichment testing, clustering, and cross-omics functional interpretation.
#'
#' @usage
#' MOLECULAR_SIGNATURES
#'
#' @format
#' An object of class \code{list}. Each element corresponds to a molecular
#' signature and contains the set of genes or features defining that pathway.
#'
#' @details
#' This object is derived from MSigDB and includes pathways from multiple
#' collections (e.g., canonical pathways, curated gene sets). No filtering or
#' pruning is applied beyond harmonization of identifiers for compatibility
#' with downstream analyses.
#'
#' @keywords datasets pathways msigdb
#'
"MOLECULAR_SIGNATURES"


#' @title Gene set ID to pathway name mapping
#'
#' @description
#' Lookup table mapping internal or database-specific gene set identifiers to
#' human-readable pathway names. This object is primarily used for labeling,
#' reporting, and cross-referencing enrichment analysis results.
#'
#' @usage
#' SET_TO_ID
#'
#' @format
#' An object of class \code{data.frame} with at least two columns:
#'
#' \describe{
#'   \item{\code{set_id}}{character; unique gene set identifier.}
#'   \item{\code{set_name}}{character; descriptive pathway or gene set name.}
#' }
#'
#' @details
#' This mapping is intended to be used alongside \code{MOLECULAR_SIGNATURES}
#' and enrichment results (e.g., CAMERA, ORA) to translate compact set IDs into
#' interpretable pathway labels for figures and tables.
#'
#' @keywords datasets annotation
#'
"SET_TO_ID"


#' @title CAMERA enrichment analysis results
#'
#' @description
#' Results from Correlation Adjusted MEan RAnk (CAMERA) enrichment analysis
#' applied to differential analysis results across tissues and molecular
#' modalities. These results are generated using the
#' \code{\link{run_cameraPR}} wrapper, which performs pre-ranked CAMERA testing
#' on z-statistics derived from differential analyses.
#'
#' CAMERA enrichment accounts for inter-gene correlation and provides a
#' competitive test of whether molecular signatures show coordinated
#' up- or down-regulation relative to the background.
#'
#' @usage
#' CAMERA_RESULTS
#'
#' @format
#' An object of class \code{data.frame} with the following columns:
#'
#' \describe{
#'   \item{\code{tissue}}{factor; the tissue in which the enrichment was tested.}
#'
#'   \item{\code{assay}}{factor; the molecular assay or omics layer.}
#'
#'   \item{\code{contrast_type}}{factor; the type of contrast (e.g., acute,
#'   chronic, sex-stratified).}
#'
#'   \item{\code{contrast}}{factor; the specific contrast tested.}
#'
#'   \item{\code{contrast_short}}{factor; abbreviated contrast label.}
#'
#'   \item{\code{collection}}{factor; the molecular signature collection.}
#'
#'   \item{\code{database}}{factor; the source database for the molecular
#'   signature (e.g., Reactome, GO, KEGG).}
#'
#'   \item{\code{set_id}}{character; unique identifier for the molecular
#'   signature.}
#'
#'   \item{\code{set}}{character; full molecular signature or pathway name.}
#'
#'   \item{\code{set_short}}{character; shortened or standardized pathway name.}
#'
#'   \item{\code{set_size}}{integer; number of molecules from the set present in
#'   the differential analysis results.}
#'
#'   \item{\code{set_size_DB}}{integer; total number of molecules in the set as
#'   defined in the source database.}
#'
#'   \item{\code{size_ratio}}{numeric; ratio of \code{set_size} to
#'   \code{set_size_DB}, providing a measure of set coverage.}
#'
#'   \item{\code{direction}}{factor; direction of enrichment, either
#'   \code{"Up"} or \code{"Down"}.}
#'
#'   \item{\code{t}}{numeric; two-sample t-statistic from CAMERA.}
#'
#'   \item{\code{df}}{integer; degrees of freedom.}
#'
#'   \item{\code{z.std}}{numeric; standard Normal equivalent of the test
#'   statistic.}
#'
#'   \item{\code{p_value}}{numeric; two-sided p-value.}
#'
#'   \item{\code{adj_p_value}}{numeric; BH-adjusted p-value, adjusted within
#'   tissue, assay, contrast, and collection.}
#' }
#'
#' @details
#' CAMERA_RESULTS summarizes pathway-level enrichment patterns across tissues,
#' assays, and experimental contrasts. These results are typically used for
#' downstream visualization (e.g., enrichment heatmaps), cluster annotation,
#' and integrative interpretation alongside
#' \code{\link{MOLECULAR_SIGNATURES}} and \code{\link{SET_TO_ID}}.
#'
#' @seealso
#' \code{\link{run_cameraPR}},
#' \code{\link{MOLECULAR_SIGNATURES}},
#' \code{\link{SET_TO_ID}}
#'
#' @keywords datasets enrichment camera

"CAMERA_RESULTS"



