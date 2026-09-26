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
  verbose = TRUE,
  load_clinical = FALSE
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
  Naming any single metabolomics platform is the same as naming
  `"metab"`: the research platforms live in one stacked object per
  tissue, so all of them are loaded and the platform is read off the
  `platform` column. `"metab-t-clinical"` is the exception — it is
  clinical chemistry, is not in the stack, and is gated by
  `load_clinical`.

- single_matrix:

  logical; if `TRUE`, returns a single combined `data.frame` across all
  selected tissues and assays. If `FALSE` (default), returns a nested
  list.

- verbose:

  logical; toggle verbosity.

- load_clinical:

  logical; whether to include the clinical chemistry omes
  ([`clinical_ome_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/clinical_ome_list.md):
  `"prot-clinical"` and `"metab-t-clinical"`). `FALSE` by default, so
  `"all"` returns the research omes and nothing changes for callers
  written before c2.0 split clinical chemistry out. Set `TRUE` to
  include them; they are dropped even when named unless it is set.

## Value

If `single_matrix = FALSE`, a nested list of `data.frame` objects. The
top-level names correspond to tissues and the second-level names to
assays — the same nesting
[`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md)
returns, so the two tiers can be walked together. The research
metabolomics platforms arrive as a single `"metab"` element per tissue
rather than one element per platform.

If `single_matrix = TRUE`, a single `data.frame` containing all selected
summary statistics, with missing columns filled as `NA`.

Each table carries `tissue`, `assay`, `randomGroupCode`, `Timepoint`,
`feature_id`, `Count`, `Mean` and `SD`. The metabolomics tables carry
`assay = "metab"` and one further column, `platform`, naming the
platform the row was measured on. That is how the `*_DA` objects are
labelled, so a join between the two tiers no longer has to translate
between two names for the same ome.

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

## Load metabolomics only. One table per tissue, every platform in it.
metab_stats = load_summary_stats(selected_omes = "metab")
unique(metab_stats[["blood"]][["metab"]][["platform"]])

## Load adipose transcriptomics as a single table
adipose_rna = load_summary_stats(
  selected_tissues = "adipose",
  selected_omes = "transcript-rna-seq",
  single_matrix = TRUE
)
} # }
```
