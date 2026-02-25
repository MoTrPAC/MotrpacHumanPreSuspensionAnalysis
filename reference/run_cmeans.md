# Fuzzy C-Means (FCM) Clustering

Fuzzy C-Means (FCM) Clustering

## Usage

``` r
run_cmeans(
  selected_tissues = c("all", "adipose", "blood", "muscle"),
  selected_omes = c("transcript-rna-seq", "prot-pr", "prot-ol", "prot-ph", "metab"),
  num_clusters_adipose = 13L,
  num_clusters_blood = 13L,
  num_clusters_muscle = 12L,
  modality = c("both", "Endur", "Resist")
)
```

## Arguments

- selected_tissues:

  character; passed to
  [`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md).
  One or more of the following: `"all"`, `"muscle"`, `"adipose"`, or
  `"blood"`.

- selected_omes:

  character; one or more of `"transcript-rna-seq"`, `"prot-pr"`,
  `"prot-ol"`, `"prot-ph"`, or `"metab"`.

- num_clusters_adipose, num_clusters_blood, num_clusters_muscle:

  integer; the number of clusters desired for each tissue. Defaults to
  13 for adipose, 13 for blood, and 12 for muscle. If more than one
  number is provided,
  [`Mfuzz::Dmin`](https://rdrr.io/pkg/Mfuzz/man/Dmin.html) will be used
  to generate a plot of the minimum centroid distances, and the user
  will be asked to specify the optimal cluster number in the console for
  each tissue.

- modality:

  character; which exercise modalities should be used for FCM? One of
  `"both"`, `"Endur"`, or `"Resist"`.

## Value

A named list of objects where names are tissues. Each object is of class
`"fclust"` with additional list component `"input"` for the matrix of
scaled z-scores used as input for FCM (both modalities included, even
when `modality != "both"`).

## See also

[`plot_cmeans`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/plot_cmeans.md),
[`run_cluster_cameraPR`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cluster_cameraPR.md),
[`run_cluster_ORA`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cluster_ORA.md)

## Author

Tyler Sagendorf, Christopher Jin

## Examples

``` r
if (FALSE) { # \dontrun{
  x1 <- run_cmeans()
  names(x1) # list available components

  # Try a range of cluster numbers for one tissue
  x2 <- run_cmeans(selected_tissues = "adipose",
                   num_clusters_adipose = 3:14)

  # FCM for a single modality
  x3 <- run_cmeans(selected_tissues = "adipose",
                   modality = "Endur")
} # }
```
