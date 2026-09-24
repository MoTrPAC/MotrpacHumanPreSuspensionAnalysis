# Transcription-Factor Phosphosite Regulator Pool

The prot-ph features whose gene is a curated human transcription factor.
The Human TFs database extract of Lambert et al (PMID:29425488) is
restricted to the entries curated as `Is TF? == "Yes"`, and those gene
symbols are mapped through
[`HUMAN_FEATURE_TO_GENE`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/HUMAN_FEATURE_TO_GENE.md)
to every prot-ph feature sharing the symbol.

This is a regulator pool, not an annotation table. It exists to restrict
a network inference to phosphosites sitting on transcription factors,
and the `feature_id` column is what such a caller subsets on. The
curated TF call, the DNA-binding domain and the rest of the Lambert
annotation are not carried here; consult the source below for those.

## Usage

``` r
UTORONTO_TFs
```

## Format

A data frame with 1,381 rows and 2 columns, one row per prot-ph feature.

- feature_id:

  character; the prot-ph phosphosite identifier, unique within the
  table.

- gene_symbol:

  character; the HGNC symbol of the transcription factor the site sits
  on. 383 distinct symbols, since a factor generally carries more than
  one site.

## Source

<https://humantfs.ccbr.utoronto.ca/download.php>

## See also

[`HUMAN_FEATURE_TO_GENE`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/HUMAN_FEATURE_TO_GENE.md)
