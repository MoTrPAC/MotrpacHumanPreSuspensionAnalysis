# Load Epigenomic Differential Analysis Results from AWS

Loads publicly released epigenomic differential analysis (DA) results
for selected tissues and omes directly from an AWS-backed content
delivery network (CDN). These datasets are hosted on Amazon CloudFront
and are accessed via HTTPS without requiring authentication or AWS
credentials. Only needed for epigenetic files because these are too
large for the package.

This function iterates over the requested tissue–ome combinations and
retrieves each corresponding DA table using
[`.load_single_ome_tissue_AWS`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/dot-load_single_ome_tissue_AWS.md).
The returned object mirrors the nested list structure used elsewhere in
the package, with tissues at the top level and omes at the second level.

The AWS-hosted files correspond to the public release of the human
pre-suspension sedentary adult epigenomics analyses.

## Usage

``` r
load_DA_from_AWS(selected_tissues, selected_omes)
```

## Arguments

- selected_tissues:

  character vector; one or more tissues for which differential analysis
  results should be retrieved.

- selected_omes:

  character vector; one or more epigenomic omes to load (e.g.,
  `"epigen-methylcap-seq"`, `"epigen-atac-seq"`).

## Value

A nested list of `data.frame` objects. The top-level names correspond to
tissues, and the second-level names correspond to omes. Each entry
contains a differential analysis table loaded from the AWS public
release.

## Details

Files are served via Amazon CloudFront and downloaded on demand at
runtime. Connection timeouts are increased internally to accommodate
large epigenomic result files.

## See also

[`.load_single_ome_tissue_AWS`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/dot-load_single_ome_tissue_AWS.md)
[`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md)

## Author

Christopher Jin

## Examples

``` r
if (FALSE) { # \dontrun{
DA_epigen <- load_DA_from_AWS(
  selected_omes = c("epigen-methylcap-seq"),
  selected_tissues = c("muscle", "adipose")
)
} # }
```
