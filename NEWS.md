
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
