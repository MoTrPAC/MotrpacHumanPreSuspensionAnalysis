# PTM-SEA Results

PTM-SEA enrichment of PTMsigDB signatures in the prot-ph
differential-analysis z-statistics (see `PTMSEA_INPUT`), run with the
Broad Institute PTM-SEA / ssGSEA2.0 container, for the acute-exercise
EE-CON and RE-CON contrasts. Muscle carries three post-exercise
timepoints crossed with the two exercise modalities; adipose was sampled
at `post_3.5_4_hr` only. Same layout as
[`CAMERA_RESULTS`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/CAMERA_RESULTS.md),
with `NES` in place of `t`, `df` and `z.std`, so it can be passed to
[`plot_enrich_heatmap`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/plot_enrich_heatmap.md).

## Usage

``` r
PTMSEA_RESULTS
```

## Format

A `data.frame` with 3,910 rows (muscle: 506 signatures x 6 contrasts;
adipose: 437 signatures x 2 contrasts) and 17 columns:

- `tissue`:

  factor; `"muscle"` or `"adipose"`.

- `assay`:

  factor; `"prot-ph"`.

- `contrast_type`:

  factor; `"exercise_with_controls"`.

- `contrast`:

  factor; the contrast, as in `CONTRAST_CONVERTER$contrast`.

- `contrast_short`:

  factor; abbreviated contrast label.

- `collection`:

  factor; `"PTMSIGDB"`.

- `database`:

  factor; the PTMsigDB signature category, the prefix of `set` (e.g.
  `"KINASE-PSP"`, `"PERT-P100-DIA2"`).

- `set_id`:

  factor; the PTMsigDB signature ID. Not an ID from
  [`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md).

- `set`:

  factor; the PTMsigDB signature ID.

- `set_short`:

  factor; signature name followed by `database` in parentheses.

- `set_size`:

  integer; number of signature sites measured.

- `set_size_DB`:

  integer; number of sites in the signature, as reported by PTM-SEA.

- `size_ratio`:

  numeric; `set_size / set_size_DB`.

- `direction`:

  factor; `"Up"` if `NES > 0`, else `"Down"`.

- `NES`:

  numeric; normalized enrichment score.

- `p_value`:

  numeric; nominal p-value.

- `adj_p_value`:

  numeric; PTM-SEA FDR, adjusted within tissue and contrast.

## Source

Muscle: precovid-analyses `figures/muscle/Figure5/` (Natalie Clark).
Adipose: Cheehoon Ahn. Provenance in `data-raw/PTMSEA/README.md`.

## See also

[`CAMERA_RESULTS`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/CAMERA_RESULTS.md),
[`plot_enrich_heatmap`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/plot_enrich_heatmap.md),
[`PTMSEA_INPUT`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/PTMSEA_INPUT.md)
