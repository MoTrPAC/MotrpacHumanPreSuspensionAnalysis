# Over-Representation Analysis (ORA)

Fast over-representation analysis.

## Usage

``` r
run_ORA(
  input,
  background,
  database = names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
  path_to_gmt = NULL,
  min_size = 5L,
  overlap_cutoff = 0.7
)
```

## Arguments

- input:

  character; vector of "interesting" features. Most likely genes, but
  may be RefMet metabolite IDs or singly-phosphorylated peptides.

- background:

  character; vector of all features that were analyzed. Used to filter
  the molecular signatures.

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

- `collection`:

  factor; the broad molecular signature collection. Only included when
  `path_to_gmt` is `NULL`. See
  [`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)
  for details.

- `database`:

  factor; the molecular signature database. Only included when
  `path_to_gmt` is `NULL`. See
  [`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)
  for details.

- `set_id`:

  character; a unique ID for the molecular signature. Only included when
  `path_to_gmt` is `NULL`. See
  [`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)
  for details.

- `set`:

  character; the molecular signature being tested. For global proteomics
  and transcriptomics, these are gene sets. For phosphoproteomics, these
  are kinase sets.

- `set_short`:

  character; a shortened version of `set`. Only included when
  `path_to_gmt` is `NULL`. See
  [`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)
  for details.

- `set_size`:

  integer; the number of molecules in the set that were present in the
  `background` vector.

- `set_size_DB`:

  integer; the number of molecules in the set, as defined in
  [`MOLECULAR_SIGNATURES`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/MOLECULAR_SIGNATURES.md).

- `size_ratio`:

  numeric; the ratio of `set_size` to `set_size_DB`, rounded to the
  nearest thousandth. A measure of confidence that the gene set being
  tested is correctly described by the entry in the `set` column. While
  smaller values do not necessarily indicate that the results are
  unreliable, terms from the gene set databases should be treated with
  caution.

- `set_size_in_input`:

  integer; the number of molecules in the set that were present in the
  `input` vector.

- `input_size`:

  integer; the size of the `input` vector.

- `background_size`:

  integer; the size of the `background` vector.

- `p_value`:

  numeric; the two-sided p-value.

- `adj_p_value`:

  numeric; the BH-adjusted p-value. P-values are adjusted within each
  combination of tissue, assay, contrast, and collection.

## See also

[`run_cluster_ORA`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cluster_ORA.md),
[`run_cameraPR`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cameraPR.md)

## Examples

``` r
# Use genes from all gene sets as the background
bg <- MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES
bg[c("PSP", "REFMET")] <- NULL
bg <- unique(unlist(bg))

# Use all genes from the MitoCarta OXPHOS term as the interesting subset
input <- MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES$MITOCARTA[["MITOCARTA_OXPHOS"]]

res <- run_ORA(input = input, background = bg)
#> Loading required package: limma
head(res)
#>   collection  database set_id
#> 1  MITOCARTA MITOCARTA  11461
#> 2         C2  REACTOME  13334
#> 3  MITOCARTA MITOCARTA  11463
#> 4         C2  REACTOME  13668
#> 5         C5      GOCC  08593
#> 6         C2        WP  14188
#>                                                                                                                         set
#> 1                                                                                                          MITOCARTA_OXPHOS
#> 2 REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS
#> 3                                                                                                 MITOCARTA_OXPHOS subunits
#> 4                                                     REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT
#> 5                                                                                             GOCC_ORGANELLE_INNER_MEMBRANE
#> 6                                                                 WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA
#>                                                      set_short set_size
#> 1                                                       OXPHOS      169
#> 2        REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP...(13334)      127
#> 3                                              OXPHOS subunits      102
#> 4 REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY...(13668)      178
#> 5                                GOCC_ORGANELLE_INNER_MEMBRANE      553
#> 6    WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA      105
#>   set_size_DB size_ratio set_size_in_input input_size background_size
#> 1         169          1               169        169           40385
#> 2         127          1               114        169           40385
#> 3         102          1               102        169           40385
#> 4         178          1               114        169           40385
#> 5         553          1               140        169           40385
#> 6         105          1                97        169           40385
#>         p_value   adj_p_value
#> 1  0.000000e+00  0.000000e+00
#> 2 5.734173e-277 9.920120e-275
#> 3 1.947040e-260 3.699375e-259
#> 4 5.434940e-245 4.701223e-243
#> 5 1.825406e-237 1.522389e-234
#> 6 3.382743e-235 1.950715e-233
```
