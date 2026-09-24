# Create a heatmap of select features or a choice pathway.

Create a heatmap for user-specified features for a given tissue and ome
combination. This function replaces `plot_pathway_features`

## Usage

``` r
plot_feature_heatmap(
  feature_ids = NULL,
  set_id = NULL,
  platforms = NULL,
  selected_tissue = c("adipose", "blood", "muscle"),
  selected_ome = c("transcript-rna-seq", "prot-pr", "prot-ph", "metab", "prot-ol"),
  contrast_type = "exercise_with_controls",
  column_title = "",
  filename,
  max_size = NULL,
  post_min = NULL,
  post_hr = NULL,
  full_modality_names = FALSE,
  multi_tissue_clust_rows = FALSE,
  right_annotation = NULL,
  heatmap_args = list(),
  draw_args = list(),
  return_drawing = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- feature_ids:

  character or `NULL`; vector of feature IDs that will appear in the
  heatmap. Must be a subset of `HUMAN_FEATURE_TO_ID[["feature_id"]]`.

- set_id:

  character of `NULL`; if not `NULL`, the feature IDs in the set will be
  used to filter rows for the heatmap, but those IDs will not be used
  for the heatmap row labels.

- platforms:

  character or `NULL`; the platform(s) used to filter metabolites. If
  `NULL` (default), metabolites in `feature_ids` will be selected from
  all available platforms. If metabolites appear in more than one
  platform, the platform will appear before the metabolite name in the
  row names of the heatmap.

- selected_tissue:

  character; the tissue that will be used to create the heatmap.

- selected_ome:

  character; the ome that will be used to create the heatmap.

- contrast_type:

  character; the type of contrasts to plot. One of
  "exercise_with_controls" (default), "exercise_no_controls",
  "Endur_vs_Resist", "baseline", or "control_only".

- column_title:

  character; the title you'd like to include for the columns. Usually
  empty

- filename:

  character; optional file name used to save the heatmap. If provided,
  the heatmap will not be drawn. Ignored when `return_drawing = TRUE`.

- max_size:

  numeric; largest number of pathways to display, if a pathway has too
  many features. Will automatically chose the first n pathways.

- post_min:

  numeric; for experiments, timepoints could be either 15, 30, or 45
  minutes depending on the analysis, this argument allows users to
  specify which of those 3 values to use, by default the value is NULL
  and will provide a generic label of post 15/30/45 min

- post_hr:

  numeric; for experiments, timepoints could be either 3.5 or 4 hours
  depending on the analysis, this argument allows users to specify which
  of those 2 values to use, by default the value is NULL and will
  provide a generic label of post 3.5/4 hr

- full_modality_names:

  logical; if TRUE the modality values are set as Endurance Exercise and
  Resistance Exercise but if FALSE the modality values are set as EE and
  RE (by default the value is FALSE)

- multi_tissue_clust_rows:

  logical; whether to cluster rows when more than one tissue is
  selected. Rows are restricted to features with a value in every
  column, since missing values break row clustering. Single-tissue
  heatmaps are always clustered.

- right_annotation:

  `NULL`, a
  [`HeatmapAnnotation`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html),
  or a function. A function receives the row labels in heatmap row order
  and must return a row `HeatmapAnnotation`; use it when the annotation
  depends on which feature each row is. A `HeatmapAnnotation` is used as
  is and must already be in row order.

- heatmap_args:

  list; arguments passed to
  [`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html). They
  override the defaults set here, e.g. `list(cluster_rows = FALSE)`.

- draw_args:

  list; arguments passed to
  [`draw`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).
  They override the defaults set here, e.g. `list(newpage = FALSE)`.

- return_drawing:

  logical; if `TRUE`, nothing is drawn or saved. Instead a list is
  returned so the caller controls the graphics device.

- verbose:

  logical; for specific warnings and additional information.

- ...:

  Additional parameters to be added to a ComplexHeatmap call

## Value

If `return_drawing = FALSE` (default), nothing; the heatmap is drawn on
the current device, or saved to `filename` if provided. If
`return_drawing = TRUE`, a list with components `draw`, a function with
no arguments that draws the heatmap on the current device without
starting a new page, and `width` and `height`, the suggested page size
in inches.

## Author

Tyler Sagendorf, Damon Leach, Christopher Jin

## Examples

``` r
if (FALSE) { # \dontrun{
plot_feature_heatmap(set_id = "11725",
  DA_list = DA_list,
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "prot-ph",
  filename = "sandbox/test_feature_heatmap.pdf")

# Draw on a device the caller opens, with a row annotation
hm <- plot_feature_heatmap(
  feature_ids = c("ENSG00000109819.9", "ENSG00000112715.26",
                  "ENSG00000119508.18", "ENSG00000162772.17"),
  selected_tissue = c("muscle", "adipose"),
  selected_ome = "transcript-rna-seq",
  multi_tissue_clust_rows = TRUE,
  right_annotation = function(row_labels) {
    ComplexHeatmap::rowAnnotation(group = rep("A", length(row_labels)))
  },
  return_drawing = TRUE)
grDevices::pdf("heatmap.pdf", width = hm$width, height = hm$height)
hm$draw()
grDevices::dev.off()
} # }
```
