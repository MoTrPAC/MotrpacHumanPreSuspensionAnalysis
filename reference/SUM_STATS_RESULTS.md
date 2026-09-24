# Summary Statistics by Tissue, Assay, and Platform

Group- and timepoint-level summary statistics computed across multiple
tissues, molecular assays, and analytical platforms. Each dataset
contains means and standard deviations calculated within randomization
groups and timepoints. These objects are intended for descriptive and
exploratory use and are not differential analysis results. Available for
all omes/tissue/platforms except methylcap, which is processed
separately.

For any datasets with missing values, the values reflect means/sds/n
with missing values excluded.

## Usage

``` r
## Adipose
ADIPOSE_EPIGEN_METHYLCAP_SEQ_SUM_STATS
ADIPOSE_METAB_SUM_STATS

ADIPOSE_PROT_PH_SUM_STATS
ADIPOSE_PROT_PR_SUM_STATS
ADIPOSE_TRANSCRIPT_RNA_SEQ_SUM_STATS

## Blood
BLOOD_EPIGEN_METHYLCAP_SEQ_SUM_STATS
BLOOD_METAB_SUM_STATS

BLOOD_PROT_OL_SUM_STATS
BLOOD_METAB_T_CLINICAL_SUM_STATS
BLOOD_PROT_CLINICAL_SUM_STATS
BLOOD_TRANSCRIPT_RNA_SEQ_SUM_STATS

## Muscle
MUSCLE_EPIGEN_METHYLCAP_SEQ_SUM_STATS
MUSCLE_METAB_SUM_STATS

MUSCLE_PROT_PH_SUM_STATS
MUSCLE_PROT_PR_SUM_STATS
MUSCLE_TRANSCRIPT_RNA_SEQ_SUM_STATS
MUSCLE_EPIGEN_ATAC_SEQ_SUM_STATS
```

## Format

Each object is a `data.frame` with one row per feature per randomization
group and timepoint, carrying the columns below in this order. They are
named and ordered to agree with the `*_DA` objects: what the rows are,
then what identifies a row, then the statistics.

- tissue:

  adipose, blood or muscle.

- assay:

  The ome. Every metabolomics platform is stacked under `"metab"` and
  named in `platform`, exactly as the `*_DA` objects do it; every other
  ome names itself here.

- platform:

  The metabolomics platform, for example `"metab-u-rppos"` or
  `"metab-t-clinical"`. Present on metabolomics objects only, which is
  also how the `*_DA` objects carry it. Objects for any other ome do not
  have this column.

- randomGroupCode:

  Randomization group: ADUControl, ADUEndur, ADUResist.

- Timepoint:

  Timepoint within the acute bout.

- feature_id:

  Feature identifier, in the same namespace as the differential analysis
  for that tissue and ome.

- Count:

  Number of samples summarised, after excluding missing values.

- Mean:

  Mean across those samples.

- SD:

  Standard deviation across those samples, `NA` when `Count` is 1.

Before v2.0.1 the platform was written into `assay`, there was no
`platform` column, and the identifying columns came first with `tissue`
and `assay` last. There was also one object per metabolomics platform
rather than one per tissue.

## Details

Datasets are stratified by tissue (e.g., adipose, blood, muscle) and
assay (metabolomics, proteomics, transcriptomics, epigenomics) — one
object per tissue and assay, named the way the `*_DA` objects are named.

Every metabolomics object reads `assay = "metab"` and names its platform
in a `platform` column (e.g. `"metab-u-rppos"`, `"metab-t-tca"`), which
is how the `*_DA` objects are labelled.

The research platforms are not one object each. They are stacked into a
single `{TISSUE}_METAB_SUM_STATS` per tissue, which is how
`{TISSUE}_METAB_DA` is keyed, so the two tiers nest the same way.
`BLOOD_METAB_T_CLINICAL_SUM_STATS` is clinical chemistry: it carries the
same labelling but stays its own object, on this tier and the
differential-analysis tier alike, because it shares five analytes with
the research platforms (Cortisol, Glucose, Glycerol, KET, NEFA) that
would otherwise sit in one object twice. `platform` is what tells the
two apart — `assay` does not.
