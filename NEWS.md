# MotrpacHumanPreSuspensionAnalysis 2.0.3

## Data objects

- `HUMAN_FEATURE_TO_GENE` gains three columns and has 13 where it had 10. No row and no
  existing column changed: all 1,920,618 rows and all ten previously shipped columns are
  identical to 2.0.2, so nothing that reads the table today reads anything different. Code
  selecting columns by position has to be updated; code selecting by name does not.

  - `custom_annotation` (factor) and `relationship_to_gene` (numeric) say where an
    ATAC-seq or MethylCap-seq peak sits relative to the gene it was assigned to — the
    region it falls in (`"Promoter (<=1kb)"`, `"Intron"`, `"Distal Intergenic"`, and seven
    others) and the signed base-pair distance to that gene, `0` where the peak overlaps it.
    Both were produced by the pipeline all along and dropped before the table was built, so
    an epigenomics feature arrived carrying only a gene: a peak in a promoter and a peak
    40 kb into an intron were indistinguishable once mapped. They are populated for all
    1,852,716 epigenomics rows and `NA` everywhere else.

    Unlike `confident_site` they are not a per-tissue measurement — they are derived from
    the peak coordinates in the `feature_id` — so they take the same value in every tissue
    a peak appears in and are not collapsed. This was checked rather than assumed: none of
    the 306,788 ATAC or 1,545,930 MethylCap feature identifiers shared between tissues
    disagree.

  - `confident_site` (logical) is the phosphosite localization flag, `NA` outside
    `prot-ph`. It is **collapsed across tissues** — `TRUE` only where a site is confidently
    localized in every tissue that measured it — because this table is keyed on
    `(assay, feature_id)` and has no tissue column, and 859 sites disagree between muscle
    and adipose. Read `*_PROT_PH_QC$feature_metadata` in
    `MotrpacHumanPreSuspensionData` when tissue-specific localization matters;
    `preprocess_PTMSEA()` already does, and is unaffected by this addition.

## Documentation

- `?HUMAN_FEATURE_TO_GENE` documented `assay` as a factor. It is a character vector, and
  has been for as long as the table has been built this way.

- The documented `assay` values did not include `"prot-clinical"`, which the 2.0 split of
  clinical chemistry into a metabolomics and a proteomics assay introduced. All eight
  values the column actually takes are now listed.

# MotrpacHumanPreSuspensionAnalysis 2.0.2

## Breaking changes

- `load_differential_analysis(epigen = TRUE)` reads the epigenomics tables from Google
  Cloud Storage via gsutil instead of the public CloudFront release, and defaults to the
  current precovid-repro staging bucket. It now requires `repo_local_dir`, and gains
  `gsutil` and `bucket` arguments; `plot_single_feature()` passes all three through.
- The unexported `load_DA_from_AWS()` and `.load_single_ome_tissue_AWS()` are removed.
  Its pinned `version = "1.2"` no longer matched the atac-seq tables, which are at v2.0.

## Other changes

- `MotrpacBicQC` returns as a `Suggests` dependency, loaded only on the `epigen = TRUE`
  path, which reports that it relies on the gsutil implementation.

# MotrpacHumanPreSuspensionAnalysis 2.0.1

## Breaking changes to data objects

- The `*_SUM_STATS` objects are named, ordered and keyed the way the `*_DA` objects are.
  Every metabolomics platform is now labelled `assay = "metab"` with the platform
  in its own `platform` column, where before the platform was written into `assay` and
  there was no `platform` column. Code that reads `assay` on a metabolomics summary
  statistic, or that joins one to a differential analysis result, has to read `platform`
  instead. Objects for every other ome are unchanged in this respect and still have no
  `platform` column.

  There are 17 objects where there were 46. The research metabolomics platforms are no
  longer one object each: they are stacked into a single `{TISSUE}_METAB_SUM_STATS` per
  tissue — `ADIPOSE_METAB_SUM_STATS` (12 platforms, 8,004 rows), `BLOOD_METAB_SUM_STATS`
  (10, 21,717), `MUSCLE_METAB_SUM_STATS` (10, 10,692) — which is how `{TISSUE}_METAB_DA`
  is keyed, so the two tiers nest the same way. Code naming one of the 32 per-platform
  objects has to read the tissue's stack and filter `platform` instead.

  `BLOOD_METAB_T_CLINICAL_SUM_STATS` is not in the stack and remains its own object. It
  is clinical chemistry, it is its own object on the differential-analysis tier too, and
  it shares five analytes with the research platforms (Cortisol, Glucose, Glycerol, KET,
  NEFA) that would otherwise sit in one object twice. It carries the same
  `assay = "metab"` / `platform = "metab-t-clinical"` labelling either way.

  The column order changed with it, from
  `randomGroupCode, feature_id, Timepoint, Count, Mean, SD, tissue, assay` to
  `tissue, assay, [platform], randomGroupCode, Timepoint, feature_id, Count, Mean, SD`.
  Code selecting columns by position has to be updated; code selecting by name does not.

  No value changed: Count, Mean and SD are identical to v2.0.0 for every row of all 46
  objects.

- `plot_single_feature()` requires these objects and errors on summary statistics that
  predate them, rather than drawing a panel with no points.

- `load_summary_stats()` returns the research metabolomics platforms as a single `"metab"`
  element per tissue rather than one element per platform — the nesting
  `load_differential_analysis()` returns, so the two tiers can be walked together. Naming
  one platform still loads them all, as before.

- `load_differential_analysis()` no longer rewrites `BLOOD_METAB_T_CLINICAL_DA`'s `assay`
  to `"metab-t-clinical"` on read. Objects are returned as stored. That rewrite existed
  because the summary statistics of the day named the platform in `assay` and the two
  tiers therefore disagreed; they now agree, so it is gone.

  The consequence is worth stating plainly: `assay` names the assay family, not the
  platform. Clinical chemistry and the research platforms both read `"metab"`, so a key
  that must separate them has to include `platform` — `(tissue, assay, feature_id)` alone
  selects two rows for each of the five shared analytes. `load_clinical = FALSE` is still
  the default, so clinical rows only arrive when asked for.

## Bug fixes

- `plot_single_feature()` draws the same legend for every tissue, whether or not that
  tissue has a timepoint below the p threshold. Combining plots with
  `patchwork::plot_layout(guides = "collect")` previously produced a repeated p
  threshold legend, because collection only merges guides that are identical and a
  tissue with nothing significant contributed a one-key legend. A single plot with no
  significant timepoints now shows both p threshold keys rather than only `adj p >=`.

- `plot_single_feature()` plots clinical chemistry only when a clinical ome is
  requested. It previously appended the clinical rows whenever the feature name
  matched an analyte, so `selected_omes = "transcript-rna-seq"` with `"Glucose"`
  returned a clinical chemistry plot. `selected_omes` now accepts the omes in
  `clinical_ome_list()` by name and `"all"` includes them, matching how
  `load_differential_analysis()` treats clinical chemistry.

  Breaking: a request that names another ome no longer returns clinical chemistry
  alongside it, and `"metab"` no longer implies `"metab-t-clinical"`. Analytes
  measured both clinically and on a research platform, such as Cortisol and
  Lactate, return only what was asked for.

- `plot_single_feature()` loads every non-epigenetic tissue and ome once and filters
  afterwards, rather than assembling the request ome by ome. Clinical chemistry is no
  longer a special case appended after the load, and the differential analysis and the
  summary statistics are put in one vocabulary before either is filtered.

  `plot_single_feature("VEGFA")` works again. The default `selected_omes = "all"` was
  broken for every non-metab feature by `filter(platform != "metab-t-conv")`: `platform`
  is NA on non-metab rows and `filter` drops NA, so the filter deleted the whole
  non-metab payload and the feature was reported as absent from the data.

- `plot_single_feature()` no longer excludes the `metab-t-conv` platform, which is now
  plotted and labelled `Conv. Metab (log2)`. It has no `assay_codes` row, so without
  that fallback its facet strip read `NA`. Note that it is the `metab-t-clinical`
  measurement on a log2 scale, so Glucose, Glycerol, KET and NEFA in blood now return a
  panel from each.


# MotrpacHumanPreSuspensionAnalysis 2.0.0

Data objects regenerated by the precovid-repro pipeline for the v2.0 release.

## New data

- BLOOD_METAB_T_CLINICAL_DA, BLOOD_METAB_T_CLINICAL_SUM_STATS, BLOOD_PROT_CLINICAL_DA, BLOOD_PROT_CLINICAL_SUM_STATS

  Clinical chemistry, one assay in v1.3, is split into a metabolomics and a proteomics assay.

## Removed data

- `BLOOD_CLINICAL_CHEMISTRY_SUM_STATS` — replaced by the v2.0 split into BLOOD_METAB_T_CLINICAL_SUM_STATS and BLOOD_PROT_CLINICAL_SUM_STATS.
- `BLOOD_EPIGEN_ATAC_SEQ_SUM_STATS` — no blood ATAC feature is significant at adj_p_value < 0.05 this cycle, and epigenomics summary statistics carry significant features only, so no object is built.
- `CLIN_CHEMISTRY_DA` — replaced by the v2.0 split into BLOOD_METAB_T_CLINICAL_DA and BLOOD_PROT_CLINICAL_DA.

## Breaking changes to data objects

- 11 objects drop the `CI.L`, `CI.R` columns; code that selects them will error: ADIPOSE_METAB_DA, ADIPOSE_PROT_PH_DA, ADIPOSE_PROT_PR_DA, ADIPOSE_TRNSCRPT_DA, BLOOD_METAB_DA, BLOOD_PROT_OL_DA, BLOOD_TRNSCRPT_DA, MUSCLE_METAB_DA, MUSCLE_PROT_PH_DA, MUSCLE_PROT_PR_DA, MUSCLE_TRNSCRPT_DA.

- HUMAN_FEATURE_TO_GENE gains the `flanking_sequence` column.

## Provenance

- `inst/PROVENANCE.tsv` records each object as regenerated, staged verbatim, or carried forward. This payload: 4 added, 2 carried forward, 8 identical, 3 removed, 12 schema change, 54 values differ.


# MotrpacHumanPreSuspensionAnalysis 0.2.4

This release addresses installation failures on R 4.6 and removes the `Remotes:`
field that was breaking dependency resolution.

Verified with clean installs on R 4.4 (Bioconductor 3.20), R 4.5 (3.22), and
R 4.6 (3.23).

# MotrpacHumanPreSuspensionAnalysis 0.2.3

## User-facing functions

- clinical chemistry is labeled `clinical-chemistry` instead of `clinical_chemistry` 
- `plot_single_feature()` now clarifies clinical chemistry features correctly; the y-axis
  label now clarifies that clinical chemistry features are shown on an absolute scale
  while all other features are shown as `log2(normalized value)`.
- Updated `HUMAN_FEATURE_TO_GENE` following changes described in 0.2.2., which now maps to the features in qc-norm properly.   

## Internals and package checks

- Namespace hygiene updates - proper imports have been labeled throughout
- Trimmed `globalVariables()`, consolidated `@importFrom` tags, added the `grid.rect`
  import, filled in the package `Description`, and added `RColorBrewer`, `doParallel`,
  `parallel`, and `randomForest` to `Suggests`. 

## Backend: data generation and release tooling (data-raw)

These changes only affect analysts with sample-level data from `MotrpacHumanPreSuspensionData`.

- Added a Google Cloud bucket validation pipeline under
  `data-raw/google_cloud_bucket_checks/`. The numbered scripts snapshot the production
  GCS bucket, copy it to an isolated staging folder, diff staging against the local
  updated files, upload changed files (removing superseded versions), and validate the
  staging structure and values against the installed package. A `run_validation_pipeline.R`
  driver orchestrates the steps; see the directory `README.md` for setup.
- Moved ATAC peak annotation into `generate_atac_qc_norm.R`
  (`.annotate_atac_features()` / `pre_cawg_get_peak_annotations_hs()`) and refined the
  gene-mapping logic. Also updated the methylcap, clinical, and metabolomics generation
  scripts.


# MotrpacHumanPreSuspensionAnalysis 0.2.2


## Backend: Feature metadata gene annotation (data-raw)

These changes will only make functional differences for analysts with sample level data from `MotrpacHumanPreSuspensionData`

- Added gene-level annotation to feature_metadata outputs for all proteomics omes. Previously,
  feature_metadata for Prot-PR, Prot-PH, and Prot-OL contained only raw provenance columns (UniProt
  accessions, PTM identifiers, redundant IDs) with no standardized gene symbol or Ensembl mapping.
  Feature metadata files now include `gene_symbol`, `ensembl_gene`, and `entrez_gene` columns
  derived from a three-round BioMart lookup strategy.
- Added `.annotate_olink()` to `generate_prot_ol_qc_norm.R`: resolves Olink protein metadata
  (UniProt accession from `uniprot_entry`) to gene symbols and Ensembl IDs. Lookup proceeds via
  UniProt → gene symbol → Entrez ID, with three progressive BioMart passes to maximize coverage.
- Added `.annotate_prot_pr()` to `generate_prot_pr_qc_norm.R`: extends the per-tissue
  `feature_metadata_output` from `PR@rdesc` with gene annotations. Uses `feature_id` (= `protein_id`,
  a UniProt accession) taken directly from the file as the BioMart lookup key.
- Added `.annotate_prot_ph()` to `generate_prot_ph_qc_norm.R`: same annotation for
  phosphoproteomics. Uses `protein_id` (not `feature_id`, which is the PTM site identifier) taken
  directly from the file as the UniProt lookup key, preserving the full PTM `feature_id` in the output.
- Prot-PR and Prot-PH annotation now resolves UniProt accessions to genes via a single BioMart
  `getBM` lookup (mirroring the Olink `.annotate_olink()` implementation), using the UniProt
  accession read directly from the source files rather than a downloaded UniProt ID mapping table.


# MotrpacHumanPreSuspensionAnalysis 0.2.1

- Added clinical chemistry differential analysis using the same structure as the molecular differential analysis;
  see: `CLIN_CHEMISTRY_DA`

# MotrpacHumanPreSuspensionAnalysis 0.2.0

## Documentation and website

- Added pkgdown GitHub Actions workflow for automated website deployment.
- Added pkgdown configuration updates for GitHub Pages publishing.
- Clarified consortium-only optional functionality in the README.

## Dependency handling

- Made `MotrpacHumanPreSuspensionData` an optional consortium-only runtime dependency
  instead of a public package dependency for CI/pkgdown resolution.

# MotrpacHumanPreSuspensionAnalysis 0.1.0

## Initial release

- Public release of summary statistics and modeling outputs from the
  MoTrPAC human pre-COVID suspension phase.
- Added functions for loading differential analysis and summary statistics.
- Added clustering, enrichment, and visualization utilities.
- Added vignettes and package website scaffolding.
