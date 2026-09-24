# List the Metabolomics Platforms

A character vector containing the targeted and untargeted metabolomics
platforms.

## Usage

``` r
metab_only_list()
```

## Value

A character vector containing the targeted and untargeted metabolomics
platforms.

## Details

These are the research metabolomics platforms.

`"metab-t-clinical"` is deliberately NOT here, though it is a targeted
metabolomics platform and a valid choice in
[`ome_available_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/ome_available_list.md).
It is clinical chemistry, and the loaders gate it behind `load_clinical`
rather than folding it into the research platforms — see
[`clinical_ome_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/clinical_ome_list.md).

## See also

[`ome_available_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/ome_available_list.md),
[`clinical_ome_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/clinical_ome_list.md)

## Examples

``` r
metab_only_list()
#>  [1] "metab-u-hilicpos"  "metab-u-ionpneg"   "metab-u-lrpneg"   
#>  [4] "metab-u-lrppos"    "metab-u-rpneg"     "metab-u-rppos"    
#>  [7] "metab-t-amines"    "metab-t-conv"      "metab-t-imm-crt"  
#> [10] "metab-t-oxylipneg" "metab-t-tca"       "metab-t-nuc"      
#> [13] "metab-t-acoa"      "metab-t-ka"       
```
