# Save a Plot of Fuzzy C-Means Cluster Trajectories

Essentially, a ggplot2 version of `Mfuzz::plot.mfuzz2` with a different
color scheme that saves a plot to a file, rather than displaying it.

## Usage

``` r
plot_cmeans(
  FCM,
  filename,
  keep_clusters,
  min_membership = 0,
  common_ylim = TRUE,
  ncol = 4L,
  dpi = 250,
  modality = c("both", "Endur", "Resist")
)
```

## Arguments

- FCM:

  one of the elements of the list produced by
  [`run_cmeans`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cmeans.md).
  See examples.

- filename:

  File name to create on disk.

- keep_clusters:

  integer; the identifier numbers of clusters to include in the plot.
  Defaults to all clusters.

- min_membership:

  numeric; to be included in the plot, the probability that a feature
  belongs to a cluster should be \\\geq\\ `min_membership`.

- common_ylim:

  logical; whether the plots for each cluster should have the same
  y-axis limits.

- ncol:

  integer; the maximum number of clusters to place on a single row. The
  number of rows is determined automatically.

- dpi:

  Plot resolution. Also accepts a string input: "retina" (320), "print"
  (300), or "screen" (72). Only applies when converting pixel units, as
  is typical for raster output types.

- modality:

  character; which exercise modalities should be used for plotting? One
  of `"both"`, `"Endur"`, or `"Resist"`.

## See also

[`run_cmeans`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cmeans.md)

## Author

Tyler Sagendorf

## Examples

``` r
if (FALSE) { # \dontrun{
# FCM results
FCM <- run_cmeans(selected_tissues = "adipose")

plot_cmeans(FCM = FCM[["adipose"]],
            filename = "adipose_clusters.pdf",
            min_membership = 0.3)
} # }
```
