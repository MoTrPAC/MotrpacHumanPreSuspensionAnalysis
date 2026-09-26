# Load Differential Analysis Results

Load Differential Analysis Results

## Usage

``` r
load_differential_analysis(
  selected_omes = "all",
  selected_tissues = "all",
  single_matrix = FALSE,
  epigen = FALSE,
  repo_local_dir = NULL,
  combine_with_featgene = FALSE,
  verbose = TRUE,
  load_clinical = FALSE
)
```

## Arguments

- selected_omes:

  character; one of
  [`ome_available_list`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/ome_available_list.md).

- selected_tissues:

  character; one of
  [`tissue_available_list`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/tissue_available_list.md).

- single_matrix:

  logical; if `TRUE`, returns a single `data.frame` containing all
  results. Otherwise, returns a list of `data.frame` objects (default).

- epigen:

  logical; a toggle of TRUE/FALSE if epigenetics data is desired. The
  epigenomics tables are not shipped in the package; they are downloaded
  from the public c2.0 release on the MoTrPAC CloudFront distribution at
  run time. No credentials are needed, but downloads are not cached and
  are slow due to file sizes.

- repo_local_dir:

  Deprecated and ignored. Epigenomics files are no longer downloaded to
  a local cache; supplying it only prints a message.

- combine_with_featgene:

  logical; whether to include columns from `HUMAN_FEATURE_TO_GENE` in
  the output.

- verbose:

  logical; whether or not to display messages for some warnings.

- load_clinical:

  logical; whether to include the clinical chemistry omes
  ([`clinical_ome_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/clinical_ome_list.md):
  `"prot-clinical"` and `"metab-t-clinical"`). `FALSE` by default, so
  `"all"` returns the research omes and nothing changes for callers
  written before c2.0 split clinical chemistry out. Set `TRUE` to
  include them; they are dropped even when named unless it is set.

## Value

A nested list of `data.table` objects. The top level names are the
tissues, while the second level names are the omes. Each table may
possess the following columns:

- tissue:

  factor; the tissue.

- assay:

  factor; the assay family, not the platform. Every metabolomics table
  reads `"metab"` — clinical chemistry included — so `assay` alone does
  not separate `metab-t-clinical` from the research platforms, and five
  analytes (Cortisol, Glucose, Glycerol, KET, NEFA) exist on both.
  Include `platform` in any key that has to tell them apart.

- platform:

  factor; (metabolomics only) the metabolomics platform.

- full_model:

  factor; full model containing predictors and any covariates.

- contrast:

  factor; full contrast (up to 33).

- contrast_short:

  factor; shortened version of the contrasts.

- contrast_type:

  factor; one of "exercise_with_controls", "exercise_no_controls",
  "Endur_vs_Resist", "baseline", or "control_only".

- contrast_category:

  factor; one of "EE-CON", "RE-CON", "EE-EE", "RE-RE", "EE-RE", or
  "CON-CON".

- feature_id:

  factor; feature identifier (may be proteins, phosphosites,
  transcripts, metabolites/lipids, peaks, or GpGs).

- logFC:

  numeric; difference between the group means in the contrast.

- CI.L_calculated:

  numeric; lower bound of the 95\\ interval on `logFC`, computed as
  `logFC - (logFC / t) * qt(0.975, df)` against that contrast's own
  residual degrees of freedom.

- CI.R_calculated:

  numeric; upper bound of the same interval. The pair is named
  `_calculated` to keep it distinct from `topTable`'s `CI.L`/`CI.R`,
  which these tables do not carry: for a `dream` fit those bound every
  contrast by the first contrast's degrees of freedom.

- degrees_of_freedom:

  numeric; the per-feature residual degrees of freedom used to shrink
  that feature's variance. Note this is NOT the degrees of freedom the
  p-value was computed from, which is the per-contrast Satterthwaite
  value plus the prior; p-values cannot be recomputed from this column.

- logLik:

  numeric; log likelihood of differential expression.

- AveExpr:

  numeric; mean of all sample-level values for that feature.

- methylation_diff:

  numeric; methylation difference.

- t:

  numeric; moderated t-statistic.

- z.std:

  numeric; standard normal equivalent of the t-statistic (z-scores).

- p_value:

  numeric; p-value.

- adj_p_value:

  numeric; p-values adjusted within each combination of tissue, assay,
  platform, and contrast using the Benjamini-Hochberg method to control
  the false discovery rate.

## Author

Tyler Sagendorf Christopher Jin

## Examples

``` r
DA_list <- load_differential_analysis() # default behavior
#> You've requested one or more epigenetic omes (via explicit selection or "all") but `epigen = FALSE`, so epigenetic data will be skipped. Set `epigen = TRUE` to load epigenetic data.
#> Clinical omes (prot-clinical, metab-t-clinical) are skipped; set `load_clinical = TRUE` to include them.
#> Please remember that the lowest CV Metabolite is chosen and the
#>             relevant refmet name is used. If you're not able to find your desired
#>             metabolite, look through the METABOLOMICS_CVS object for the relevant
#>             refmet/feature name.

# Structure of a single object
str(DA_list[["adipose"]][["prot-pr"]])
#> Classes ‘data.table’ and 'data.frame':   63504 obs. of  20 variables:
#>  $ tissue            : chr  "adipose" "adipose" "adipose" "adipose" ...
#>  $ assay             : chr  "prot-pr" "prot-pr" "prot-pr" "prot-pr" ...
#>  $ full_model        : Factor w/ 1 level "~ 0 + group_timepoint +  BMI + calculatedAge + Sex + (1 | pid)": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ contrast          : Factor w/ 9 levels "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3"| __truncated__,..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ contrast_short    : Factor w/ 33 levels "Endur.during_20_min - Control.during_20_min (delta-delta)",..: 5 5 5 5 5 5 5 5 5 5 ...
#>  $ contrast_type     : Factor w/ 5 levels "exercise_with_controls",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ contrast_category : Factor w/ 6 levels "EE-CON","RE-CON",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ randomGroupCode   : chr  "ADUEndur" "ADUEndur" "ADUEndur" "ADUEndur" ...
#>  $ Timepoint         : Factor w/ 7 levels "pre_exercise",..: 6 6 6 6 6 6 6 6 6 6 ...
#>  $ feature_id        : chr  "O00287" "Q8TE02" "Q8NHG8" "P10912" ...
#>  $ logFC             : num  0.74 0.578 0.609 -0.857 -0.71 ...
#>  $ CI.L_calculated   : num  0.533 0.414 0.463 -1.087 -0.933 ...
#>  $ CI.R_calculated   : num  0.946 0.742 0.756 -0.627 -0.487 ...
#>  $ degrees_of_freedom: num  27 24 16.1 22.7 24 ...
#>  $ logLik            : num  11.02 15.87 20.89 2.42 6.32 ...
#>  $ t                 : num  7.31 7.21 8.75 -7.73 -6.52 ...
#>  $ AveExpr           : num  0.0188 -1.1554 -0.2189 -0.205 -0.5905 ...
#>  $ z.std             : num  5.53 5.38 5.35 -5.3 -5.04 ...
#>  $ p_value           : num  3.22e-08 7.53e-08 8.96e-08 1.14e-07 4.62e-07 ...
#>  $ adj_p_value       : num  0.000202 0.000202 0.000202 0.000202 0.000651 ...
#>  - attr(*, ".internal.selfref")=<pointer: (nil)> 
#>  - attr(*, "sorted")= chr [1:4] "full_model" "contrast" "p_value" "feature_id"

# Un-nest list
DA_list <- unlist(DA_list, recursive = FALSE)
names(DA_list)
#>  [1] "adipose.metab"              "adipose.prot-ph"           
#>  [3] "adipose.prot-pr"            "adipose.transcript-rna-seq"
#>  [5] "blood.metab"                "blood.prot-ol"             
#>  [7] "blood.transcript-rna-seq"   "muscle.metab"              
#>  [9] "muscle.prot-ph"             "muscle.prot-pr"            
#> [11] "muscle.transcript-rna-seq" 

# Include epigen data
if (FALSE) { # \dontrun{
DA_list <- load_differential_analysis(epigen = TRUE)
} # }
```
