# Gene set ID to pathway name mapping

Lookup table mapping internal or database-specific gene set identifiers
to human-readable pathway names. This object is primarily used for
labeling, reporting, and cross-referencing enrichment analysis results.

## Usage

``` r
SET_TO_ID
```

## Format

An object of class `data.frame` with at least two columns:

- `set_id`:

  character; unique gene set identifier.

- `set_name`:

  character; descriptive pathway or gene set name.

## Details

This mapping is intended to be used alongside `MOLECULAR_SIGNATURES` and
enrichment results (e.g., CAMERA, ORA) to translate compact set IDs into
interpretable pathway labels for figures and tables.
