# CAMERA enrichment analysis results

Results from Correlation Adjusted MEan RAnk (CAMERA) enrichment analysis
applied to differential analysis results across tissues and molecular
modalities. These results are generated using the
[`run_cameraPR`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cameraPR.md)
wrapper, which performs pre-ranked CAMERA testing on z-statistics
derived from differential analyses.

CAMERA enrichment accounts for inter-gene correlation and provides a
competitive test of whether molecular signatures show coordinated up- or
down-regulation relative to the background.

## Usage

``` r
CAMERA_RESULTS
```

## Format

An object of class `data.frame` with the following columns:

- `tissue`:

  factor; the tissue in which the enrichment was tested.

- `assay`:

  factor; the molecular assay or omics layer.

- `contrast_type`:

  factor; the type of contrast (e.g., acute, chronic, sex-stratified).

- `contrast`:

  factor; the specific contrast tested.

- `contrast_short`:

  factor; abbreviated contrast label.

- `collection`:

  factor; the molecular signature collection.

- `database`:

  factor; the source database for the molecular signature (e.g.,
  Reactome, GO, KEGG).

- `set_id`:

  character; unique identifier for the molecular signature.

- `set`:

  character; full molecular signature or pathway name.

- `set_short`:

  character; shortened or standardized pathway name.

- `set_size`:

  integer; number of molecules from the set present in the differential
  analysis results.

- `set_size_DB`:

  integer; total number of molecules in the set as defined in the source
  database.

- `size_ratio`:

  numeric; ratio of `set_size` to `set_size_DB`, providing a measure of
  set coverage.

- `direction`:

  factor; direction of enrichment, either `"Up"` or `"Down"`.

- `t`:

  numeric; two-sample t-statistic from CAMERA.

- `df`:

  integer; degrees of freedom.

- `z.std`:

  numeric; standard Normal equivalent of the test statistic.

- `p_value`:

  numeric; two-sided p-value.

- `adj_p_value`:

  numeric; BH-adjusted p-value, adjusted within tissue, assay, contrast,
  and collection.

## Details

CAMERA_RESULTS summarizes pathway-level enrichment patterns across
tissues, assays, and experimental contrasts. These results are typically
used for downstream visualization (e.g., enrichment heatmaps), cluster
annotation, and integrative interpretation alongside
[`MOLECULAR_SIGNATURES`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/MOLECULAR_SIGNATURES.md)
and
[`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md).

## See also

[`run_cameraPR`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/run_cameraPR.md),
[`MOLECULAR_SIGNATURES`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/MOLECULAR_SIGNATURES.md),
[`SET_TO_ID`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/SET_TO_ID.md)
