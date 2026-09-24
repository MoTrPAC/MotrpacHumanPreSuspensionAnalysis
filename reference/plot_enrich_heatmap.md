# Create a bubble heatmap of molecular signatures from CAMERA-PR or PTM-SEA results

Create a heatmap of molecular signatures for a given tissue, ome, and
contrast type combination from CAMERA-PR results (`CAMERA_RESULTS` or
the output of
[`run_cameraPR()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cameraPR.md))
or PTM-SEA results (`PTMSEA_RESULTS`).

## Usage

``` r
plot_enrich_heatmap(
  x,
  set_ids = NULL,
  selected_tissues = c("adipose", "blood", "muscle"),
  selected_ome = c("transcript-rna-seq", "prot-pr", "prot-ph", "metab", "prot-ol"),
  contrast_type = c("exercise_with_controls", "exercise_no_controls", "Endur_vs_Resist",
    "baseline", "control_only"),
  padj_cutoff = 0.05,
  n_top = 6L,
  zscore_colors = c("#3366ff", "darkred"),
  filename = NULL,
  return_drawing = FALSE
)
```

## Arguments

- x:

  a `data.frame` of CAMERA-PR results with a `z.std` column, such as
  `CAMERA_RESULTS` or the output of
  [`run_cameraPR()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cameraPR.md),
  or PTM-SEA results with an `NES` column, such as `PTMSEA_RESULTS`.

- set_ids:

  character or `NULL`; one or more of `x$set_id`. The dataset `x` will
  be filtered to these sets. If `NULL`, the top 6 most significant sets
  from each contrast will be selected.

- selected_tissues:

  character; one or more tissues that will be used to create the
  heatmap. If multiple tissues are selected, only those terms tested in
  at least two tissues will be plotted.

- selected_ome:

  character; the ome that will be used to create the heatmap.

- contrast_type:

  character; the type of contrasts to plot. One of
  "exercise_with_controls" (default), "exercise_no_controls",
  "Endur_vs_Resist", "baseline", or "control_only".

- padj_cutoff:

  numeric; the p-value cutoff to add an asterisk indicating statistical
  significance.

- n_top:

  integer; if `set_ids` is `NULL`, the `n_top` most significant
  molecular signatures from each tissue and contrast will be selected.

- zscore_colors:

  character; length 2 vector of colors for the smallest and largest
  z-scores. Default is "#3366ff" (blue) and "darkred".

- filename:

  character; path to a PDF file to save the heatmap. Ignored if
  `return_drawing = TRUE`.

- return_drawing:

  logical; if `TRUE`, nothing is drawn or saved. Instead a list is
  returned so the caller controls the graphics device.

## Value

If `return_drawing = FALSE` (default), nothing; the heatmap is saved to
`filename`. If `return_drawing = TRUE`, a list with components `draw`, a
function with no arguments that draws the heatmap on the current device
without starting a new page; `width` and `height`, the suggested page
size in inches; and `n_sets`, the number of sets drawn.

## Author

Tyler Sagendorf
