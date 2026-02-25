# Changelog

## MotrpacHumanPreSuspensionAnalysis 0.2.0

### Documentation and website

- Added pkgdown GitHub Actions workflow for automated website
  deployment.
- Added pkgdown configuration updates for GitHub Pages publishing.
- Clarified consortium-only optional functionality in the README.

### Dependency handling

- Made `MotrpacHumanPreSuspensionData` an optional consortium-only
  runtime dependency instead of a public package dependency for
  CI/pkgdown resolution.

## MotrpacHumanPreSuspensionAnalysis 0.1.0

### Initial release

- Public release of summary statistics and modeling outputs from the
  MoTrPAC human pre-COVID suspension phase.
- Added functions for loading differential analysis and summary
  statistics.
- Added clustering, enrichment, and visualization utilities.
- Added vignettes and package website scaffolding.
