# Load Summary Statistics for Normalized Expression Data

Loads group- and timepoint-level summary statistics for normalized
expression data across tissues, molecular assays, and analytical
platforms. Summary statistics consist of means and standard deviations
computed within randomization groups and timepoints. Limited to the
initial acute bout, training data is excluded here.

These datasets are intended for descriptive and exploratory analyses.
Sample-level data are not distributed with this package and are
available upon request via <https://motrpac-data.org>.

For metabolomics assays, summary statistics are computed after filtering
redundant metabolites; see the Methods section of the associated
documentation for details.

## Usage

``` r
load_summary_stats(
  selected_tissues = "all",
  selected_omes = "all",
  single_matrix = FALSE,
  verbose = TRUE
)
```

## Arguments

- selected_tissues:

  character; tissues to include. One or more of `"adipose"`, `"blood"`,
  `"muscle"`, or `"all"`.

- selected_omes:

  character; molecular assays to include. One or more of
  `"transcript-rna-seq"`, `"prot-pr"`, `"prot-ph"`, `"prot-ol"`,
  `"epigen-atac-seq"`, `"epigen-methylcap-seq"`, `"metab"`, or `"all"`.
  Selecting `"metab"` loads all metabolomics platforms.

- single_matrix:

  logical; if `TRUE`, returns a single combined `data.frame` across all
  selected tissues and assays. If `FALSE` (default), returns a nested
  list.

- verbose:

  logical; toggle verbosity.

## Value

If `single_matrix = FALSE`, a nested list of `data.frame` objects. The
top-level names correspond to tissues, and the second-level names
correspond to assays or platforms.

If `single_matrix = TRUE`, a single `data.frame` containing all selected
summary statistics, with missing columns filled as `NA`.

## Details

Summary statistics were filtered to only those that qualified for
differential analysis. This means for proteomics/phosphoproteomics,
samples required a paired n\>=3 to be included. See the methods in the
manuscript for more information. In epigenetic assays, features were
filtered for significant features (FDR\<0.05) for file size purposes.

## Author

Christopher Jin

## Examples

``` r
if (FALSE) { # \dontrun{
## Load all summary statistics
sum_stats = load_summary_stats()

## Load metabolomics only
metab_stats = load_summary_stats(selected_omes = "metab")

## Load adipose transcriptomics as a single table
adipose_rna = load_summary_stats(
  selected_tissues = "adipose",
  selected_omes = "transcript-rna-seq",
  single_matrix = TRUE
)
} # }
```
