# Wrapper for CAMERA-PR

A wrapper for
[`cameraPR.matrix`](https://rdrr.io/pkg/TMSig/man/cameraPR.matrix.html)
that performs pre-ranked Correlation Adjusted MEan RAnk molecular
signature analysis of differential analysis results.

## Usage

``` r
run_cameraPR(
  DA_list = NULL,
  selected_omes = c("transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol", "metab"),
  selected_tissues = "all",
  database = setdiff(names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
    "PTMSIGDB"),
  path_to_gmt = NULL,
  min_size = 5L,
  overlap_cutoff = 0.7
)
```

## Arguments

- DA_list:

  list; a named list of `data.frame` objects, each containing
  differential analysis results for a specific tissue/assay combination.
  The list must be nested, with tissues at the top level and omes within
  tissues. If `NULL` (default), the differential analysis results will
  be generated with
  [`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md).
  Unless wishing to analyze DA results that are not in
  MotrpacHumanPreSuspensionAnalysis, this should remain `NULL`.

- selected_omes:

  character; one or more character strings selected from the following
  options: `"transcript-rna-seq"`, `"prot-pr"`, `"prot-ph"`, and
  `"metab"` (all metabolomics platforms). Passed to
  [`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md).

- selected_tissues:

  character; passed to
  [`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md).
  One or more of the following: `"all"`, `"muscle"`, `"adipose"`, or
  `"blood"`.

- database:

  character; one or more names specifying the database(s) to test.
  Options are (case insensitive) `"BIOCARTA"`, `"KEGG_MEDICUS"`,
  `"PID"`, `"REACTOME"`, `"WP"` (WikiPathways database), `"GOBP"`,
  `"GOCC"`, `"GOMF"`, `"MITOCARTA"` (MitoCarta3.0 database), `"PSP"`
  (PhosphoSitePlus kinases; only valid when `selected_omes` contains
  `"prot-ph"`), or `"REFMET"` (RefMet chemical subclasses; only valid
  when `selected_omes` contains `"metab"`). See
  [`MOLECULAR_SIGNATURES`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/MOLECULAR_SIGNATURES.md)
  for details.

- path_to_gmt:

  character; (optional) path to one or more GMT files. Passed to
  [`TMSig::readGMT`](https://rdrr.io/pkg/TMSig/man/readGMT.html). If
  provided, `database` is ignored.

- min_size:

  integer; the minimum set size for testing.

- overlap_cutoff:

  numeric; the minimum proportion of genes in each set that must appear
  in a given dataset. Used to pre-filter sets. Does not affect `"metab"`
  or `"prot-ph"` results. This will always be 0.1 for `"prot-ol"`
  results.

## Value

An object of class `data.frame` with the following columns:

`tissue`

:   factor; the tissue.

## See also

[`MOLECULAR_SIGNATURES`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/MOLECULAR_SIGNATURES.md),
[`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)

## Author

Tyler Sagendorf

## Examples

``` r
if (FALSE) { # \dontrun{
# Test Reactome and Gene Ontology Biological Processes
# databases on muscle global proteomics.
x1 <- run_cameraPR(selected_omes = "prot-pr",
                   selected_tissues = "muscle",
                   database = c("GOBP", "REACTOME"))

# Test all omes using all molecular signature databases
x2 <- run_cameraPR()
} # }
```
