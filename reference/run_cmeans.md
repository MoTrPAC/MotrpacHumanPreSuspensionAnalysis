# Fuzzy C-Means (FCM) Clustering

Fuzzy c-means clustering of features by the shape of their z-score
trajectories across the `"exercise_with_controls"` contrasts (the two
`during` contrasts excluded). The pre-computed
[`FCM_CLUSTERS`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/FCM_RESULTS.md)
was built with the default arguments.

## Usage

``` r
run_cmeans(
  DA_list = NULL,
  selected_tissues = c("all", "adipose", "blood", "muscle"),
  selected_omes = c("transcript-rna-seq", "prot-pr", "prot-ol", "prot-ph", "metab"),
  num_clusters_adipose = 13L,
  num_clusters_blood = 12L,
  num_clusters_muscle = 12L,
  modality = c("both", "Endur", "Resist")
)
```

## Arguments

- DA_list:

  list; a named list of `data.frame` objects, each containing
  differential analysis results for a specific tissue/assay combination.
  The list is nested, with tissues at the top level and omes within
  tissues, or already flattened with names of the form `"tissue.assay"`.
  If `NULL` (default), the differential analysis results will be
  generated with
  [`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md).
  Unless wishing to analyze DA results that are not in
  MotrpacHumanPreSuspensionAnalysis, this should remain `NULL`.

- selected_tissues:

  character; passed to
  [`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md).
  One or more of the following: `"all"`, `"muscle"`, `"adipose"`, or
  `"blood"`.

- selected_omes:

  character; one or more of `"transcript-rna-seq"`, `"prot-pr"`,
  `"prot-ol"`, `"prot-ph"`, or `"metab"`.

- num_clusters_adipose, num_clusters_blood, num_clusters_muscle:

  integer; the number of clusters for each tissue, a single value each.
  Defaults to 13 for adipose, 12 for blood, and 12 for muscle, the
  values used for `FCM_CLUSTERS`. They were chosen from a sweep of
  cluster numbers in the motrpac-human-presuspension-repro pipeline
  (step 13), which plots the minimum centroid distance and related
  diagnostics for each tissue.

- modality:

  character; which exercise modalities should be used for FCM? One of
  `"both"`, `"Endur"`, or `"Resist"`.

## Value

A named list of objects where names are tissues. Each object is of class
`"fclust"` with additional list component `"input"` for the matrix of
scaled z-scores used as input for FCM (both modalities included, even
when `modality != "both"`).

## See also

[`FCM_CLUSTERS`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/FCM_RESULTS.md),
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

  # Reuse differential analysis results already in memory
  DA_list <- load_differential_analysis(selected_tissues = "adipose")
  x2 <- run_cmeans(DA_list = DA_list,
                   selected_tissues = "adipose")

  # FCM for a single modality
  x3 <- run_cmeans(selected_tissues = "adipose",
                   modality = "Endur")
} # }
```
