# List the Available Tissues

Vector of tissues available in the analysis. These are the "full names"
for each of the tissues.

## Usage

``` r
tissue_available_list(verbose = TRUE)
```

## Arguments

- verbose:

  logical; whether to output messages.

## Value

A character vector of the available tissues: "adipose", "blood", and
"muscle".

## Examples

``` r
tissue_available_list()
#> The available tissues are placed into overarching categories. For example, an assay using PBMCs would be categorized as blood.
#> [1] "adipose" "blood"   "muscle" 
```
