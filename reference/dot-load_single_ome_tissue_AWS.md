# Load a Single Epigenomic DA Table from AWS

Downloads and loads a single epigenomic differential analysis (DA)
result table for a specified tissue and ome from the public AWS release.
Files are hosted on Amazon CloudFront and accessed via a static HTTPS
URL.

This function is intended for internal use and is called by
[`load_DA_from_AWS`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_DA_from_AWS.md)
to populate tissue–ome combinations.

## Usage

``` r
.load_single_ome_tissue_AWS(
  selected_ome,
  selected_tissue,
  data_category = "da",
  version = "1.2"
)
```

## Arguments

- selected_ome:

  character; epigenomic ome to load (e.g., `"epigen-methylcap-seq"` or
  `"epigen-atac-seq"`).

- selected_tissue:

  character; tissue identifier.

- data_category:

  character; analysis category, default is `"da"` for differential
  analysis.

- version:

  character; dataset version identifier appended to the file name

## Value

A `data.frame` containing the differential analysis results for the
specified tissue and ome. If no matching tissue–ome combination is
found, the function returns `NULL`.

## Details

The underlying files are part of the public AWS release of the human
pre-suspension sedentary adult epigenomics analyses. File naming
conventions vary slightly across epigenomic assays and are handled
internally.

To accommodate large epigenomic result files, the global R connection
timeout is temporarily increased within this function.

## Author

Christopher Jin
