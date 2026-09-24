# List All Omes

Return a character vector containing all available omes. These are the
"full names" for each of the assays.

## Usage

``` r
ome_available_list()
```

## Value

A character vector containing all available omes.

## Details

This vector is the vocabulary every accessor gates on: `load_qc()`,
[`load_summary_stats()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_summary_stats.md)
and
[`load_differential_analysis()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md)
all resolve `"all"` through it and reject anything absent from it, so an
ome missing here is an ome no accessor can return. It must therefore
track what the pipeline actually builds, which is recorded in
`OME_TISSUE_CODE`.

Two clinical omes are included as of v2.0. Clinical chemistry was a
single `"clinical-chemistry"` assay in v1.3 and is now split into
`"metab-t-clinical"` and `"prot-clinical"`, each with its own QC,
differential-analysis and summary-statistic objects.

`"metab-meta-reg"` was removed. It was dropped from the pipeline
entirely, so no object was ever produced under that name and every
accessor offered it as a choice that returned nothing.

The `lab-*` tiers in `OME_TISSUE_CODE` are deliberately absent. They are
the raw clinical laboratory inputs, distributed as the `cln_chemistry_*`
tables rather than as omes.

## See also

[`metab_only_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/metab_only_list.md),
[`tissue_available_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/tissue_available_list.md)

## Examples

``` r
ome_available_list()
#>  [1] "prot-ol"              "prot-ph"              "prot-pr"             
#>  [4] "prot-clinical"        "transcript-rna-seq"   "epigen-methylcap-seq"
#>  [7] "epigen-atac-seq"      "metab-u-hilicpos"     "metab-u-ionpneg"     
#> [10] "metab-u-lrpneg"       "metab-u-lrppos"       "metab-u-rpneg"       
#> [13] "metab-u-rppos"        "metab-t-amines"       "metab-t-conv"        
#> [16] "metab-t-imm-crt"      "metab-t-oxylipneg"    "metab-t-tca"         
#> [19] "metab-t-nuc"          "metab-t-acoa"         "metab-t-ka"          
#> [22] "metab-t-clinical"    
```
