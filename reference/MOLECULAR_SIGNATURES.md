# Molecular pathway signatures

Collection of molecular pathway and gene set signatures derived from
MSigDB. This dataset aggregates curated pathway definitions spanning
multiple biological databases and is used as the reference universe for
gene set–based analyses throughout the package.

Molecular signatures are provided in a standardized format to support
enrichment testing, clustering, and cross-omics functional
interpretation.

## Usage

``` r
MOLECULAR_SIGNATURES
```

## Format

An object of class `list`. Each element corresponds to a molecular
signature and contains the set of genes or features defining that
pathway.

## Details

This object is derived from MSigDB and includes pathways from multiple
collections (e.g., canonical pathways, curated gene sets). No filtering
or pruning is applied beyond harmonization of identifiers for
compatibility with downstream analyses.
