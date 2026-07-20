# Assay code lookup table (vendored copy of MotrpacBicQC::assay_codes)

Lookup table mapping MoTrPAC assay codes to their human-readable names,
abbreviations, and display colors. It is used to convert internal assay
identifiers into presentation-ready labels in figures.

## Usage

``` r
assay_codes
```

## Format

A `data.frame` with 44 rows and 11 columns:

- `omics_code`:

  character; molecular modality code (e.g. `"transcriptomics"`,
  `"epigenomics"`).

- `submission_code`:

  character; code used at data submission (e.g. `"rna-seq"`,
  `"methylcap-seq"`).

- `assay`:

  character; assay identifier (e.g. `"transcript-rna-seq"`).

- `assay_code`:

  character; assay code used for joins against the `assay` column of
  differential analysis and summary stat tables.

- `assay_name`:

  character; assay name as reported by the assay site (e.g. `"RRBS"`,
  `"Methyl Capture"`).

- `cas_code`:

  character; chemical analysis site responsible for the assay (e.g.
  `"stanford"`, `"broad_prot"`).

- `assay_abbreviation`:

  character; short uppercase abbreviation (e.g. `"TRNSCRPT"`, `"ATAC"`).

- `assay_hex_colour`:

  character; hexadecimal color code used to encode the assay in figures.

- `assay_short_text`:

  character; short display label (e.g. `"RNA-seq"`, `"RRBS"`).

- `omics_text`:

  character; display label for the molecular modality (e.g.
  `"Transcriptomics"`).

- `ome_text`:

  character; display label for the ome (e.g. `"Transcriptome"`,
  `"Epigenome"`).

## Source

`MotrpacBicQC::assay_codes`, via `data-raw/assay_codes.R`.

## Details

This object is a verbatim copy of `MotrpacBicQC::assay_codes`, taken as
a snapshot by `data-raw/assay_codes.R`. It is vendored here as a
temporary measure, and
[`plot_single_feature`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/plot_single_feature.md)
now reads this copy rather than the upstream object.

The motivation for vendoring is that MotrpacBicQC imports `inspectdf`,
which has been archived on CRAN. That made MotrpacBicQC – and therefore
this package – uninstallable from a current CRAN snapshot without
pulling `inspectdf` from the CRAN Archive. Vendoring this table removed
the last use of MotrpacBicQC on an exported code path, allowing it to be
dropped from `Imports`.

Because this is a snapshot rather than a live reference, it can drift
from upstream: any assay added or relabelled in MotrpacBicQC will not
appear here until the copy is refreshed. Re-run `data-raw/assay_codes.R`
to update it. This vendoring is intended to be reverted in favor of
depending on MotrpacBicQC directly once the upstream `inspectdf`
dependency is resolved.

## See also

[`plot_single_feature`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/plot_single_feature.md),
[`OME_TISSUE_CODE`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/OME_TISSUE_CODE.md)
