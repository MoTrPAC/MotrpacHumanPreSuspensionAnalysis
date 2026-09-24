# PTM-SEA Input Matrices

Phosphoproteomics differential-analysis results formatted as the input
to PTM-SEA: one row per confidently localized phosphosite, one column
per acute-exercise contrast, values are the differential-analysis
z-statistics. This is an input file, not an enrichment result. PTM-SEA
itself is run from the Broad Institute PTM-SEA / ssGSEA2.0 container
against PTMsigDB, which is outside this package.

Sites are restricted to those flagged `confident_site` in the prot-ph
feature metadata (see `ADIPOSE_PROT_PH_QC`, `MUSCLE_PROT_PH_QC`).
Contrasts are `contrast_type == "exercise_with_controls"` with
`contrast_category` in `"EE-CON"` or `"RE-CON"`. Muscle carries three
post-exercise timepoints crossed with the two exercise modalities;
adipose was sampled at `post_3.5_4_hr` only.

Each element is a `GCT` object, the S4 class defined by the `cmapR`
package. `cmapR` is in Suggests, not Imports: the object loads without
it and its slots are reachable as `@mat`, `@rdesc`, `@cdesc`, `@rid` and
`@cid`, but anything that dispatches on the class — printing the object,
the `cmapR` accessors, `write_gct()` — needs `cmapR` installed, as does
`R CMD check` to inspect the data.

## Usage

``` r
PTMSEA_INPUT
```

## Format

A named list of two `cmapR` `GCT` objects, `"muscle"` (16,469 sites x 6
contrasts) and `"adipose"` (7,287 sites x 2 contrasts). The `mat` slot
holds the z-statistics. The `rdesc` slot holds the prot-ph feature
annotation and the per-contrast `logFC`, `p_value` and `adj_p_value`.

## Source

<https://motrpac-data.org/>
