
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
  a UniProt accession) as the lookup key.
- Added `.annotate_prot_ph()` to `generate_prot_ph_qc_norm.R`: same three-round annotation for
  phosphoproteomics. Uses `protein_id` (not `feature_id`, which is the PTM site identifier) as the
  UniProt lookup key, preserving the full PTM `feature_id` in the output.
- Added `.get_uniprot_mapping()` helper to `generate_prot_ol_qc_norm.R`: downloads the UniProt
  human ID mapping table to `tempdir()` (cached within the session) and returns the four columns
  needed for annotation (`UniProtKB-AC`, `Ensembl`, `UniProtKB-ID`, `GeneID (EntrezGene)`).


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
