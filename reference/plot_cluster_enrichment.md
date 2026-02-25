# Bubble heatmap of CAMERA-PR or ORA results from fuzzy c-means clustering

Generate a bubble heatmap from the results of nonparametric CAMERA-PR or
ORA applied to the output of fuzzy c-means clustering.

## Usage

``` r
plot_cluster_enrichment(
  x,
  set_ids = NULL,
  selected_tissues = c("adipose", "blood", "muscle"),
  selected_ome = c("transcript-rna-seq", "prot-pr", "prot-ph", "metab"),
  padj_cutoff = 0.05,
  n_top = 6L,
  column_title = NULL,
  filename
)
```

## Arguments

- x:

  a `data.frame`. Either
  [`FCM_CAMERA`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/FCM_RESULTS.md)
  or
  [`FCM_ORA`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/FCM_RESULTS.md).

- set_ids:

  character or `NULL`; one or more of `x$set_id`. The dataset `x` will
  be filtered to these sets. If `NULL`, the `n_top` most significant
  sets from each contrast will be selected.

- selected_tissues:

  character; one or more tissues that will be used to create the
  heatmap. If multiple tissues are selected, only those terms tested in
  at least two tissues will be plotted.

- selected_ome:

  character; the ome that will be used to create the heatmap.

- padj_cutoff:

  numeric; cutoff in terms of p-value for an asterisk to be included

- n_top:

  integer; if `set_ids` is `NULL`, the `n_top` most significant
  molecular signatures from each tissue and cluster will be selected.

- column_title:

  character or `NULL`; title for the heatmap.

- filename:

  character; path to a PDF file to save the heatmap.

## Value

Nothing. The heatmap is saved to a PDF specified by `filename`.

## Author

Tyler Sagendorf
