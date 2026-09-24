# Human Feature to Gene Conversion Table

A `data.table` used to convert from feature IDs to various identifiers.

## Usage

``` r
HUMAN_FEATURE_TO_GENE
```

## Format

A sorted `data.table` with 1,920,618 rows and 12 columns:

- assay:

  character; the assay (ome). One of "epigen-atac-seq",
  "epigen-methylcap-seq", "metab", "prot-clinical", "prot-ol",
  "prot-ph", "prot-pr", or "transcript-rna-seq".

- feature_id:

  factor; the feature identifier.

- entrez_gene:

  factor; Entrez gene identifier.

- gene_symbol:

  factor; gene symbol.

- ensembl_gene:

  factor; ensembl gene identifier.

- uniprot:

  factor; UniProt identifier.

- refmet_name:

  factor; RefMet metabolite name, from the pinned RefMet snapshot.

- refmet_id:

  factor; RefMet metabolite identifier, from the pinned RefMet snapshot.

- kegg_id:

  factor; KEGG identifier, from the pinned RefMet/KEGG snapshot.

- custom_annotation:

  factor; the region of the assigned gene the peak falls in — one of
  "Promoter (\<=1kb)", "Promoter (1-2kb)", "5' UTR", "3' UTR", "Exon",
  "Intron", "Overlaps Gene", "Upstream (\<5kb)", "Downstream (\<5kb)" or
  "Distal Intergenic". `NA` unless assay is "epigen-atac-seq" or
  "epigen-methylcap-seq".

- relationship_to_gene:

  numeric; the signed distance in base pairs from the peak to the
  assigned gene, `0` where the peak overlaps it. `NA` unless assay is
  "epigen-atac-seq" or "epigen-methylcap-seq".

- flanking_sequence:

  factor; flanking sequence. Only applicable if assay is "prot-ph".

## Source

Built by Stage 1 step 07 of the motrpac-human-presuspension-repro
pipeline. The `refmet_name`, `refmet_id` and `kegg_id` columns come from
a pinned offline RefMet/KEGG snapshot rather than a live Metabolomics
Workbench query, so the mapping does not move with those databases. The
peak annotations in `custom_annotation` and `relationship_to_gene` are
produced by `ChIPseeker` against the pinned Ensembl v105 `TxDb`.

## Details

`custom_annotation` and `relationship_to_gene` describe where an
ATAC-seq or MethylCap-seq peak sits relative to the gene it was assigned
to. Without them a peak in a promoter and a peak 40 kb into an intron
are indistinguishable once mapped, since both carry only the gene. Both
are derived from the peak coordinates in the `feature_id` rather than
measured per tissue, so they take the same value in every tissue a peak
appears in.

Per-tissue measurements are deliberately not carried here. This table is
keyed on `(assay, feature_id)` with no tissue column, so it could only
hold a collapse across tissues. The phosphosite localization flag
(`confident_site`) is the case in point: it was added in 2.0.3 and
removed in 2.0.7 because muscle and adipose disagree on 859 of their
7,865 shared prot-ph sites. Read it per tissue from
`*_PROT_PH_QC$feature_metadata` in MotrpacHumanPreSuspensionData.
