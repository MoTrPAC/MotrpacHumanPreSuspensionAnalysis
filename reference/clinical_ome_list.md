# List the Clinical Omes

The omes carrying clinical chemistry rather than a research assay. Every
loader takes a `load_clinical` argument that gates these, and it is
`FALSE` by default.

## Usage

``` r
clinical_ome_list()
```

## Value

A character vector of the clinical omes.

## Details

v1.3 carried clinical chemistry as a single `"clinical-chemistry"`
assay; v2.0 splits it into a metabolomics and a proteomics assay, each
with its own QC, differential-analysis and summary-statistic objects.

They are gated rather than simply included because they are a different
kind of measurement from the research omes, and every caller written
before v2.0 that asks for `"all"` is summarising the molecular
landscape. Adding them by default changes those results silently: the
clinical metabolomics differential-analysis rows overlap the combined
`*_METAB_DA` table on five analytes (Cortisol, Glycerol, KET, NEFA and
Glucose), and a caller that maps an assay to a display name puts
clinical chemistry into a row that was never meant to hold it.

Pass `load_clinical = TRUE` to get them.

## See also

[`ome_available_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/ome_available_list.md),
[`metab_only_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/metab_only_list.md)

## Examples

``` r
clinical_ome_list()
#> [1] "prot-clinical"    "metab-t-clinical"
```
