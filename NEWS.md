# MotrpacHumanPreSuspensionAnalysis 0.2.2

## Backend: QC normalization pipeline updates (data-raw)

- Added `qc_norm_visualization_helpers.R` and `qc_norm_visualization.Rmd` to compare outputs from each
  `generate_*_qc_norm()` function against the corresponding objects in `MotrpacHumanPreSuspensionData`,
  covering qc_norm feature/sample overlap, value correlation, feature_metadata ID overlap, and
  sample_metadata column-level value diffs.
- Prot-PR/PH: updated covariates file so `tmt_plex` is used as the primary covariate throughout,
  replacing the previous `plex_site` variable (an interaction of `Plex` and `Cas` specific to muscle).
  No change to differential analysis results; standardizes covariate representation across tissues.
- Prot-PR/PH: renamed feature identifier column from `protein_id` to `feature_id` in feature_metadata
  output, matching the qc_norm matrix and all other omes. feature_metadata now only includes features
  retained in the final qc_norm matrix.
- Prot-PR/PH: updated `generate_prot_pr_qc_norm()` and `generate_prot_ph_qc_norm()` to reference the
  `study` column in eQC metadata (previously `protocol`), reflecting an upstream rename in the eQC file.

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
