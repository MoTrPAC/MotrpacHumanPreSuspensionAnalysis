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
  the heatmap will not be drawn.

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

  logical; if you'd like to cluster rows specifically in a multi-tissue
  case

- verbose:

  logical; for specific warnings and additional information.

- ...:

  Additional parameters to be added to a ComplexHeatmap call

## Value

Nothing. A heatmap is drawn or saved to a file if `filename` is
provided.

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
} # }
```
