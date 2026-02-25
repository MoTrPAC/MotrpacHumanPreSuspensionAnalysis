# Analyze Fuzzy C-Means Clustering Results with CAMERA-PR

Test for localization of molecular signatures to cluster centroids.
Performs modified upper-tail signed rank tests with
[`TMSig::cameraPR.matrix`](https://rdrr.io/pkg/TMSig/man/cameraPR.matrix.html)
on the matrices of cluster membership probabilities obtained from fuzzy
c-means (FCM) clustering. Identifies molecular signatures that closely
follow the trajectories of each cluster's centroid.

## Usage

``` r
run_cluster_cameraPR(
  FCM,
  selected_omes = c("transcript-rna-seq", "prot-pr", "prot-ph", "metab"),
  selected_tissues = "all",
  database = names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
  path_to_gmt = NULL,
  min_size = 5L,
  overlap_cutoff = 0.7
)
```

## Arguments

- FCM:

  fuzzy c-means (FCM) clustering results. Output of
  [`run_cmeans`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cmeans.md).
  This should be a named list of objects of class `fclust`.

- selected_omes:

  character; one or more character strings selected from the following
  options: `"transcript-rna-seq"`, `"prot-pr"`, `"prot-ph"`, and
  `"metab"` (all metabolomics platforms). Passed to
  [`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md).

- selected_tissues:

  character; passed to
  [`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md).
  One or more of the following: `"all"`, `"muscle"`, `"adipose"`, or
  `"blood"`.

- database:

  character; one or more names specifying the database(s) to test.
  Options are (case insensitive) `"BIOCARTA"`, `"KEGG_MEDICUS"`,
  `"PID"`, `"REACTOME"`, `"WP"` (WikiPathways database), `"GOBP"`,
  `"GOCC"`, `"GOMF"`, `"MITOCARTA"` (MitoCarta3.0 database), `"PSP"`
  (PhosphoSitePlus kinases; only valid when `selected_omes` contains
  `"prot-ph"`), or `"REFMET"` (RefMet chemical subclasses; only valid
  when `selected_omes` contains `"metab"`). See
  [`MOLECULAR_SIGNATURES`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/MOLECULAR_SIGNATURES.md)
  for details.

- path_to_gmt:

  character; (optional) path to one or more GMT files. Passed to
  [`TMSig::readGMT`](https://rdrr.io/pkg/TMSig/man/readGMT.html). If
  provided, `database` is ignored.

- min_size:

  integer; the minimum set size for testing.

- overlap_cutoff:

  numeric; the minimum proportion of genes in each set that must appear
  in a given dataset. Used to pre-filter sets. Does not affect `"metab"`
  or `"prot-ph"` results. This will always be 0.1 for `"prot-ol"`
  results.

## Value

An object of class `data.frame` with the following columns:

- `tissue`:

  factor; the tissue.

- `assay`:

  factor; the omics assay.

- `cluster`:

  factor; the cluster number. The clusters have the same meaning across
  omics assays measured in the same tissue, but they have no such
  relationship across tissues.

- `collection`:

  factor; the broad molecular signature collection. See
  [`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)
  for details.

- `database`:

  factor; the molecular signature database. See
  [`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)
  for details.

- `set_id`:

  character; a unique ID for the molecular signature. See
  [`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)
  for details.

- `set`:

  character; the molecular signature being tested. For global proteomics
  and transcriptomics, these are gene sets. For phosphoproteomics, these
  are kinase sets.

- `set_short`:

  character; a shortened version of `set`. See
  [`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)
  for details.

- `set_size`:

  integer; the number of molecules in the set that were present in the
  DA results for that specific tissue/assay combination.

- `set_size_DB`:

  integer; the number of molecules in the set, as defined in the GMT
  file.

- `size_ratio`:

  numeric; the ratio of `set_size` to `set_size_DB`, rounded to the
  nearest thousandth. A measure of confidence that the gene set being
  tested is correctly described by the entry in the `set` column. While
  smaller values do not necessarily indicate that the results are
  unreliable, terms from the gene set databases should be treated with
  caution.

- `p_value`:

  numeric; the upper-tail p-value.

- `adj_p_value`:

  numeric; the BH-adjusted p-value. P-values are adjusted within each
  combination of tissue, assay, collection, and cluster.

## See also

[`run_cmeans`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cmeans.md),
[`run_cluster_ORA`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cluster_ORA.md),
[`MOLECULAR_SIGNATURES`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/MOLECULAR_SIGNATURES.md),
[`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)

## Author

Tyler Sagendorf

## Examples

``` r
if (FALSE) { # \dontrun{
  FCM <- run_cmeans()

  # Run CAMERA-PR with all available molecular signatures
  cluster_res <- run_cluster_cameraPR(FCM = FCM)
  head(cluster_res)
} # }
```
