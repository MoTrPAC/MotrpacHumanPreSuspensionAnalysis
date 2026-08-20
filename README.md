
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MotrpacHumanPreSuspensionAnalysis

<!-- badges: start -->
<!-- badges: end -->

## Overview

This R package provides the public release of summary statistics,
differential analysis results, and downstream modeling outputs from the
**Molecular Transducers of Physical Activity Consortium (MoTrPAC)**
human pre-COVID suspension cohort.

The first human cohort of MoTrPAC enrolled sedentary adults prior to
study suspension during the COVID-19 pandemic (N=175), randomized to
endurance exercise (EE), resistance exercise (RE), or non-exercise
control (CON). **This package focuses on the acute exercise bout** from
that cohort.

Participants were randomized in an approximate 8:8:3 ratio to EE, RE, or
CON groups and also to temporal profiles of biospecimen collection. A
non-exercising group was deemed critical to control for the molecular
effects of circadian rhythm, fasting, tissue sampling, and any other
non-exercise intervention stimulus. The majority of participants were
female (72%). Mean age was 41 ± 15 years, the average BMI was 26.9 ± 4.0
kg/m², and average VO2peak was 24 ± 7.0 ml/kg/min. See the MoTrPAC
manuscripts for full cohort details.

There is a larger cohort of subjects being analyzed by the MoTrPAC
Consortium for recruitment following the COVID suspension, and that
analysis will cover many more details about subgroup differences,
including information about response to longitudinal training,
heterogeneity, etc.

### What is included

- **Differential analysis** results (all tissues and omic platforms)
- **Group-level summary statistics** (n, mean, SD per
  feature/group/tissue)
- **Enrichment results** (CAMERA-PR pathway analysis)
- **Fuzzy c-means clustering** and cluster-level enrichment
- **Feature-to-gene mapping** (Ensembl v105 / GENCODE 39)
- Visualization functions for heatmaps, PCA, enrichment, and
  single-feature plots
- **Splicing analysis** significant results (FDR \< 0.05)

### What is NOT included

To protect participant privacy and comply with data-use governance
policies, individual-level (subject-level) molecular or phenotypic data
are **not** included. Such data are available only through formal data
access requests to the MoTrPAC consortium.

### Versioning note

The functions in version 0.2.0 were those used to generate the initial
bioRxiv pre-print. Over the course of reviews and additional analysis,
modifications may occur—refer to previous release history if you need to
exactly recreate pre-print Figures. We will aim to provide version
information via the GitHub “Releases” section for major version
milestones.

## How this repo fits with the others

The Pre-Suspension human work is split across four repositories. Individual-level molecular
and phenotypic data cannot be distributed publicly, so they live in a separate, access-gated
package, and everything that *can* be released publicly (aggregate results, all analysis
code) is kept clear of them.

```
                MoTrPAC BIC — consortium GCS buckets
                     (raw assay data, gated)
                                │
                         precovid-repro
        normalizes omics data, applies statistical models,
        builds every data object, versions it, uploads it,
                and carries it into both packages
                                │
              ┌─────────────────┴─────────────────┐
              ▼                                   ▼
 MotrpacHumanPreSuspensionData      MotrpacHumanPreSuspensionAnalysis
 subject-level data — access-gated  aggregate results — public
              └─────────────────┬─────────────────┘
                                ▼
               motrpac-human-presuspension-acute
               manuscript figure code + QC vignettes
                                ▼
                          manuscripts
```

| Repository | What it holds | Access |
|---|---|---|
| `precovid-repro` | the end-to-end rebuild pipeline and its pinned software environment | private; a full run needs consortium bucket access |
| [`MotrpacHumanPreSuspensionData`](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionData) | subject-level molecular and phenotypic data objects | formal data-access request to the consortium |
| [`MotrpacHumanPreSuspensionAnalysis`](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis) | differential analysis, group summary statistics, enrichment, clustering, feature-to-gene map, plotting functions | public |
| [`motrpac-human-presuspension-acute`](https://github.com/MoTrPAC/motrpac-human-presuspension-acute) | per-manuscript figure code and QC vignettes | code public; some panels need Data access |

------------------------------------------------------------------------

# Installation

## Requirements

**R \>= 4.4 is required.** Several dependencies (e.g. TMSig) require R
4.4 or later. The package will not install on older R versions.

## Important: Bioconductor Dependencies

This package relies on several Bioconductor packages
(e.g. ComplexHeatmap, Mfuzz, Biobase, TMSig). You must set the correct
Bioconductor version for your R installation before installing. Using
the wrong Bioconductor version will cause dependency failures.

| R version | Bioconductor version |
|-----------|----------------------|
| R 4.4.x   | 3.20                 |
| R 4.5.x   | 3.22                 |
| R 4.6.x   | 3.23                 |

**For R 4.4:**

``` r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.20")
if (!require("pak", quietly = TRUE)) install.packages("pak")
pak::pak("MoTrPAC/MotrpacHumanPreSuspensionAnalysis")
```

**For R 4.5:**

``` r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.22")
if (!require("pak", quietly = TRUE)) install.packages("pak")
pak::pak("MoTrPAC/MotrpacHumanPreSuspensionAnalysis")
```

**For R 4.6:**

``` r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.23")
if (!require("pak", quietly = TRUE)) install.packages("pak")
pak::pak("MoTrPAC/MotrpacHumanPreSuspensionAnalysis")
```

You can check your R version with `R.version.string` and your
Bioconductor version with `BiocManager::version()`.

We recommend building the vignettes when installing (use
`vignette(package = "MotrpacHumanPreSuspensionAnalysis")` to browse
them).

The vignette `package_overview` describes how to use the package
functions in detail, but the most common functions for those just
looking to view the results are also described here in this README.

## Troubleshooting

- **macOS (especially Apple Silicon):** If you see compilation errors or
  missing tools, install Xcode Command Line Tools first:

  ``` bash
  xcode-select --install
  ```

- **Bioconductor errors** (e.g. version mismatch warnings or failed
  binaries): Re-run `BiocManager::install(version = "X.YZ")` with the
  correct version for your R installation (see table above) before the
  GitHub install.

After installation, load the package:

``` r
library(MotrpacHumanPreSuspensionAnalysis)
```

## Consortium-only optional features

Most functionality in this package is fully public and works without
additional access. Some advanced workflows rely on the private package
`MotrpacHumanPreSuspensionData` that is available only to MoTrPAC
consortium members.

At the moment, the primary functions with this optional dependency are:

- `run_SCION()`
- `plot_precovid_cca()`

If the private package is not installed, these functions will return a
clear error message with access guidance.

## Getting help

For questions, bug reporting, and data requests for this package, please
<a
href="https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/issues"
target="_blank">submit a new issue</a> and include as many details as
possible.

------------------------------------------------------------------------

# Usage

## Omic modeling summary statistics (Differential Analysis)

See more documentation via `?load_differential_analysis` or the vignette

Note that the public release of epigenetic files is through AWS CDN, and
the file sizes are quite a bit larger than the other omes, so setting
`epigen = TRUE` can be very slow!

``` r
differential_analysis = load_differential_analysis(
  selected_omes = "all",
  selected_tissues = "all",
  single_matrix = FALSE,
  epigen = FALSE,
  combine_with_featgene = FALSE,
  verbose = TRUE
)
#> Please remember that the lowest CV Metabolite is chosen and the
#>             relevant refmet name is used. If you're not able to find your desired
#>             metabolite, look through the METABOLOMICS_CV object for the relevant
#>             refmet/feature name.
names(differential_analysis)
#> [1] "adipose" "blood"   "muscle"
names(differential_analysis[["blood"]])
#> [1] "metab"              "prot-ol"            "transcript-rna-seq"
```

By default,`load_differential_analysis` loads in the dataset in a nested
list first by tissue, then by ome. Choose whichever tissues or omes
you’d like via `selected_omes` or `selected_tissues`. You can find
available tissues via `tissue_available_list()` or
`ome_available_list()`. Or if you enter in a wrong mistaken tissue/ome,
a warning or error will help.

If you would instead like to stack the matrixes more easily, use the
`single_matrix` function, which basically unlists the list and sticks
everything into a data.frame object.

``` r
single_matrix = load_differential_analysis(single_matrix = TRUE)
#> Please remember that the lowest CV Metabolite is chosen and the
#>             relevant refmet name is used. If you're not able to find your desired
#>             metabolite, look through the METABOLOMICS_CV object for the relevant
#>             refmet/feature name.
colnames(single_matrix)
#>  [1] "tissue"             "assay"              "platform"          
#>  [4] "full_model"         "contrast"           "contrast_short"    
#>  [7] "contrast_type"      "contrast_category"  "randomGroupCode"   
#> [10] "Timepoint"          "feature_id"         "logFC"             
#> [13] "CI.L"               "CI.R"               "degrees_of_freedom"
#> [16] "logLik"             "t"                  "AveExpr"           
#> [19] "z.std"              "p_value"            "adj_p_value"
```

For a quick explanation of each of the columns, you can find this via
`?load_differential_analysis`

Importantly, this loads the differential analysis for each of
comparisons mentioned in the methods, including the comparison between
the endurance or resistance group relative to time, fasting, biopsy,
etc. matched controls, comparison between the endurance and resistance
groups directly, and finally comparisons within group without a matched
control.

The majority of the analysis is done via exercise groups relative to the
controls (“exercise_with_controls”). Make sure you filter to whichever
category you prefer before continuing with analysis.

``` r
single_matrix %>% dplyr::pull(contrast_type) %>% unique()
#> [1] exercise_with_controls exercise_no_controls   Endur_vs_Resist       
#> [4] baseline               control_only          
#> 5 Levels: exercise_with_controls exercise_no_controls ... control_only
```

If you’d like to display things in terms of the specific groups being
compared instead, you can use the ‘contrast_category’ column. (EE-CON,
RE-CON would be subsets of the ‘exercise_with_controls’ category from
above, for example.)

``` r
single_matrix %>% dplyr::pull(contrast_category) %>% unique()
#> [1] EE-CON  RE-CON  EE-EE   RE-RE   EE-RE   CON-CON
#> Levels: EE-CON RE-CON EE-EE RE-RE EE-RE CON-CON
```

The splicing data was processed in a separate analysis effort, but
significant results (FDR \< 0.05) are available in this R package as
well. The full set is available on the motrpac data hub, but is not
included here because of file size limitations. See “Exercise modulation
of the alternative splicing landscape in human tissues” for more
information.

``` r
names(SPLICING_DA)
#> [1] "adipose" "blood"   "muscle"
head(SPLICING_DA$adipose$`AS-rMATS`, 3)
#>                                           feature  Estimate Std. Error      df
#>                                            <char>     <num>      <num>   <num>
#> 1:     SE:9:36376127-36390467:36390616-36424613:-  1.766129  0.3706165 69.0000
#> 2: SE:8:105634357-105788718:105788924-105798724:+ -1.632226  0.3146590 69.0000
#> 3:         SE:7:7567484-7567607:7567688-7572436:+ -6.588839  0.1474998 24.9892
#>       t value     Pr(>|t|) log2_FoldChange   diff_psi  tissue randomGroupCode
#>         <num>        <num>           <num>      <num>  <char>          <char>
#> 1:   4.765382 1.009691e-05      0.29934027  0.1482379 adipose        ADUEndur
#> 2:  -5.187284 2.035829e-06     -0.08428790 -0.0567500 adipose        ADUEndur
#> 3: -44.670163 2.319668e-25     -0.05350315 -0.0243000 adipose        ADUEndur
#>    timepoint_baseline  timepoint_select AS_type      p_value  adj_p_value
#>                <char>            <char>  <char>        <num>        <num>
#> 1:       pre_exercise post_15_30_45_min      SE 1.009691e-05 2.341846e-02
#> 2:       pre_exercise post_15_30_45_min      SE 2.035829e-06 5.437266e-03
#> 3:       pre_exercise post_15_30_45_min      SE 2.319668e-25 6.814876e-21
#>            gene_id    assay                                           contrast
#>             <char>   <char>                                             <char>
#> 1: ENSG00000137075 AS-rMATS ADUEndur.post_15_30_45_min - ADUEndur.pre_exercise
#> 2: ENSG00000169946 AS-rMATS ADUEndur.post_15_30_45_min - ADUEndur.pre_exercise
#> 3: ENSG00000164654 AS-rMATS ADUEndur.post_15_30_45_min - ADUEndur.pre_exercise
```

## Omic modeling summary statistics

See more documentation via `?load_summary_stats`

``` r
summary_stats = load_summary_stats(
  selected_omes = "all",
  selected_tissues = "all",
  single_matrix = FALSE,
  verbose = TRUE
)
#> Only features qualifying for diffential analysis are included. For proteomics and phosphoproteomics, this means some samples with missingness patterns that lead to paired n < 3 for any group are not included here.
#> Epigenetics summary stats are trimmed to only show significant features due to file size limitations
names(summary_stats)
#> [1] "adipose" "blood"   "muscle"
names(summary_stats[["blood"]])
#>  [1] "epigen-atac-seq"      "epigen-methylcap-seq" "metab-t-amines"      
#>  [4] "metab-t-conv"         "metab-t-oxylipneg"    "metab-t-tca"         
#>  [7] "metab-u-hilicpos"     "metab-u-ionpneg"      "metab-u-lrpneg"      
#> [10] "metab-u-lrppos"       "metab-u-rpneg"        "metab-u-rppos"       
#> [13] "prot-ol"              "transcript-rna-seq"
```

By default, load_summary_stats() loads group- and timepoint-level
summary statistics for normalized expression data in a nested list
structure, organized identically to the differential-analysis datasets.
The top level corresponds to tissues, and the second level corresponds
to molecular assays or platforms.

You may subset the data using selected_tissues and selected_omes.
Available options can be queried via `tissue_available_list()` and
`ome_available_list()`. If an invalid tissue or assay is supplied,
informative warnings or errors are raised to guide correction.

If a stacked representation is preferred, set `single_matrix = TRUE`.
This unlists the nested structure and returns a single `data.frame` with
all selected tissues and assays.

Again - note that the sample level data is available for researchers
upon request via the Motrpac Consortium.

``` r
single_matrix = load_summary_stats(single_matrix = TRUE)
#> Only features qualifying for diffential analysis are included. For proteomics and phosphoproteomics, this means some samples with missingness patterns that lead to paired n < 3 for any group are not included here.
#> Epigenetics summary stats are trimmed to only show significant features due to file size limitations
colnames(single_matrix)
#> [1] "randomGroupCode" "feature_id"      "Timepoint"       "Count"          
#> [5] "Mean"            "SD"              "tissue"          "assay"
```

Summary statistics were filtered to only those that qualified for
differential analysis. This means for proteomics/phosphoproteomics,
samples required a paired n\>=3 to be included. See the methods in the
manuscript for more information.

For metabolomics assays, summary statistics are computed after filtering
redundant metabolites. Details of this filtering procedure are described
in the Methods section of the manuscript.

For epigenetic assays (ATAC, methyl), only significant features are
included, due to file size limitations.

## Enrichment Results

``` r
colnames(CAMERA_RESULTS)
#>  [1] "tissue"         "assay"          "contrast_type"  "contrast"      
#>  [5] "contrast_short" "collection"     "database"       "set_id"        
#>  [9] "set"            "set_short"      "set_size"       "set_size_DB"   
#> [13] "size_ratio"     "direction"      "t"              "df"            
#> [17] "z.std"          "p_value"        "adj_p_value"
```

Quick summary: CAMERA-PR is a method of enrichment that incorporates all
features to generate a comparison of the test statistics between
in-pathway vs out-of-pathway test statistics to see if the statistics
within pathway are significant.

This file structure is more or less just an enrichment level match for
the comparisons described in the single-matrix differential analysis
results, where all tissues and assays are included in all the analysis.

## Feature to gene file

``` r
head(HUMAN_FEATURE_TO_GENE)
#> Key: <assay, feature_id>
#>              assay               feature_id entrez_gene gene_symbol
#>             <fctr>                   <fctr>      <fctr>      <fctr>
#> 1: epigen-atac-seq chr1:100006105-100007013       23443     SLC35A3
#> 2: epigen-atac-seq chr1:100009408-100009608       23443     SLC35A3
#> 3: epigen-atac-seq   chr1:10001014-10001214      116362        RBP7
#> 4: epigen-atac-seq chr1:100010489-100010728       23443     SLC35A3
#> 5: epigen-atac-seq chr1:100021498-100021698       23443     SLC35A3
#> 6: epigen-atac-seq chr1:100024572-100024772       23443     SLC35A3
#>       ensembl_gene custom_annotation relationship_to_gene uniprot refmet_name
#>             <fctr>            <fctr>                <num>  <fctr>      <fctr>
#> 1: ENSG00000117620            Intron                    0    <NA>        <NA>
#> 2: ENSG00000117620              Exon                    0    <NA>        <NA>
#> 3: ENSG00000162444            Intron                    0    <NA>        <NA>
#> 4: ENSG00000117620            Intron                    0    <NA>        <NA>
#> 5: ENSG00000117620            Intron                    0    <NA>        <NA>
#> 6: ENSG00000117620            3' UTR                    0    <NA>        <NA>
#>    kegg_id flanking_sequence
#>     <fctr>            <fctr>
#> 1:    <NA>              <NA>
#> 2:    <NA>              <NA>
#> 3:    <NA>              <NA>
#> 4:    <NA>              <NA>
#> 5:    <NA>              <NA>
#> 6:    <NA>              <NA>
```

The feature-to-gene map links each feature tested in differential
analysis to a gene, using Ensembl version 105 (mapped to GENCODE 39) as
the gene identifier source. Proteomics feature IDs (UniProt IDs) were
mapped to gene symbols and Entrez IDs using UniProt’s mapping files.
Epigenomics features were mapped to the nearest gene using the
ChIPseeker::annotatePeak() function with Homo sapiens Ensembl release
105 gene annotations. Gene symbols, Entrez IDs, and Ensembl IDs were
assigned to features using biomaRt version 2.58.2 (Bioconductor 3.18).
This file links all of the features included in any ome/tissue in our
analysis. Use this to see how some levels of omic analysis (e.g. ATAC,
RNAseq) may link up in terms of ome names.

# For Developers

## Testing installation with Docker

A `Dockerfile` is included to verify that the package and all its
dependencies install correctly in a clean Linux environment. This is
useful for catching missing system libraries or silent dependency
failures before release.

``` bash
# Build the image (installs all Imports, Suggests, and Remotes)
docker build -t motrpac-presuspension-test .

# Run the container to confirm the package loads
docker run --rm motrpac-presuspension-test
```

To render vignettes inside the container:

``` bash
docker run --rm -v "$(pwd)/vignette_output:/output" motrpac-presuspension-test bash -c "
  apt-get update && apt-get install -y --no-install-recommends pandoc &&
  R -e \"
    rmarkdown::render('vignettes/package_overview.Rmd', output_dir='/output');
    rmarkdown::render('vignettes/differential_analysis.Rmd', output_dir='/output')
  \"
"
```

## Acknowledgements

MoTrPAC is supported by the National Institutes of Health (NIH) Common
Fund through cooperative agreements managed by the National Institute of
Diabetes and Digestive and Kidney Diseases (NIDDK), National Institute
of Arthritis and Musculoskeletal Diseases (NIAMS), and National
Institute on Aging (NIA).

Specifically, the MoTrPAC Study is supported by NIH grants U24OD026629
(Bioinformatics Center), U24DK112349, U24DK112342, U24DK112340,
U24DK112341, U24DK112326, U24DK112331, U24DK112348 (Chemical Analysis
Sites), U01AR071133, U01AR071130, U01AR071124, U01AR071128, U01AR071150,
U01AR071160, U01AR071158 (Clinical Centers), U24AR071113 (Consortium
Coordinating Center), U01AG055133, U01AG055137 and U01AG055135
(PASS/Animal Sites).

## Data Use Agreement

Recipients and their Agents agree that in publications using **any**
data from MoTrPAC public-use data sets they will acknowledge MoTrPAC as
the source of data, including the version number of the data sets used,
e.g.:

- Data used in the preparation of this article were obtained from the
  Molecular Transducers of Physical Activity Consortium (MoTrPAC)
  database, which is available for public access at
  [motrpac-data.org](motrpac-data.org).
- Data used in the preparation of this article were obtained from the
  Molecular Transducers of Physical Activity Consortium (MoTrPAC)
  Pre-CovidSuspension Data release version 1.3.0.
