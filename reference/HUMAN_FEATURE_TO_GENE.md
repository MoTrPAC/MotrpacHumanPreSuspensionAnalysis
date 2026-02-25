# Human Feature to Gene Conversion Table

A `data.table` used to convert from feature IDs to various identifiers.

## Usage

``` r
HUMAN_FEATURE_TO_GENE
```

## Format

A sorted `data.table` with 1,958,568 rows and 11 columns:

- assay:

  factor; the assay (ome). One of "epigen-atac-seq",
  "epigen-methylcap-seq", "metab", "prot-ol", "prot-ph", "prot-pr", or
  "transcript-rna-seq".

- feature_id:

  factor; the feature identifier.

- entrez_gene:

  factor; Entrez gene identifier.

- gene_symbol:

  factor; gene symbol.

- ensembl_gene:

  factor; ensembl gene identifier.

- custom_annotation:

  factor; custom feature annotation. Only applicable if assay is
  "epigen-atac-seq" or "epigen-methylcap-seq".

- relationship_to_gene:

  numeric; only applicable if assay is "epigen-atac-seq" or
  "epigen-methylcap-seq".

- uniprot:

  factor; UniProt identifier.

- refmet_name:

  factor; RefMet metabolite identifier.

- kegg_id:

  factor; KEGG identifier.

- flanking_sequence:

  factor; flanking sequence. Only applicable if assay is "prot-ph".
