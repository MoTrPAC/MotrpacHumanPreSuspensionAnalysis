# Google Cloud Bucket Validation Pipeline

This directory contains scripts for validating, staging, and deploying updated data files to
the MoTrPAC human pre-COVID suspension Google Cloud Storage (GCS) buckets. 

Note that this is currently implemented just for qc-norm and relevant sample and feature metadata. For any expected changes in the differential analysis, that will come later. 

The pipeline bridges upstream data updates (from the BIC data hub) to downstream package releases in both
the `MotrpacHumanPreSuspensionData` (data package) and `MotrpacHumanPreSuspensionAnalysis`
(analysis package).

---

## Overview

Upstream data files are versioned under the production bucket:

```
gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v{version}/
```

Before any new version of either package is released, the pipeline:

1. Snapshots the current state of the production bucket
2. Copies it to an isolated staging bucket folder
3. Diffs the staging state against a local reference folder of updated files
4. Uploads changed files to staging and removes superseded versions
5. Validates the staging folder structure and runs value checks against the installed package

The scripts in this folder implement each step. They should be run in order, or
orchestrated from the top-level `run_validation_pipeline.R` driver (see below).

---

## Prerequisites

- `gsutil` on PATH (authenticated via `gcloud auth application-default login`)
- R packages: `dplyr`, `magrittr`, `tools`, `base64enc`, `MotrpacHumanPreSuspensionData`
- Both downstream package repos cloned locally; paths set in `~/config.json`:
  - `precovid_repo_path` — path to the precovid-analyses repo (used to derive `LOCAL_UPDATED_DIR`)
  - `data_package_path` — path to `MotrpacHumanPreSuspensionData`
  - (analysis package path is resolved automatically via `here::here()`)

---

## Configuration

All bucket paths and version strings are centralized in `config.R` (sourced by every
step script) rather than hard-coded:

```r
# config.R — edit before each release cycle
PRODUCTION_BUCKET = "gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3"
STAGING_BUCKET    = "gs://pre-cawg/staging_20260511"
CURRENT_VERSION   = "1.3"
NEW_VERSION       = "1.4"

# ~/config.json keys used:
#   precovid_repo_path  — precovid-analyses repo root
#   data_package_path   — path to MotrpacHumanPreSuspensionData repo
# ANALYSIS_PKG_REPO_PATH is resolved automatically via here::here()
LOCAL_UPDATED_DIR = file.path(PRECOVID_REPO_PATH, "data", "tmp", "freeze")
```

Output directories (`snapshots/`, `diffs/`, `logs/`) are created automatically by `config.R`.

---

## Step-by-Step

### Step 1 — Snapshot the production bucket (`01_snapshot_production.R`)

**Goal:** Capture a complete, timestamped inventory of every file currently in the
production bucket so that every downstream comparison has a stable baseline.

**What it does:**

- Runs `gsutil ls -R {PRODUCTION_BUCKET}` and captures stdout
- Filters to actual file paths (lines ending in `.txt`, `.csv`, `.txt.gz`, `.rda`, `.html`)
- For each file, records GCS path, file size, and MD5 hash (from `gsutil stat`)
- Saves the result to `snapshots/snapshot_{timestamp}.tsv` with columns:
  `gcs_path | size_bytes | md5 | snapshot_timestamp`
- Copies the latest snapshot to `snapshots/latest.tsv` for use by downstream steps

**Why this matters:** Without a frozen baseline, diff comparisons between production and
local updates are unreliable if anything changes in the bucket mid-run.

---

### Step 2 — Copy production to staging (`02_copy_to_staging.R`)

**Goal:** Create an isolated working copy of the production bucket in the staging folder
so that all subsequent changes can be tested without touching the live production data.

**What it does:**

- Runs `gsutil rsync -r {PRODUCTION_BUCKET}/ {STAGING_BUCKET}/` (single-threaded; `-m` is
  intentionally omitted to avoid rate-limiting failures on large GCS-to-GCS transfers)
- Verifies the staging file count is at least as large as the snapshot from Step 1 —
  aborts if staging has fewer files than expected
- Writes a one-line confirmation (`staging_copy_verified: TRUE/FALSE`) to `logs/copy_log_*.txt`

**Why this matters:** All structural validation and upload testing happens against staging,
keeping production intact until the full pipeline passes.

---

### Step 3 — Diff staging against local updated files (`03_diff_local_vs_staging.R`)

**Goal:** Identify exactly which files have changed or been added between the incoming
local updates (`LOCAL_UPDATED_DIR`) and what currently lives in staging, with explicit
before/after version tracking for every matched file.

**What it does:**

Files are excluded from the diff in two phases before any join occurs:

**Phase 1 — exact-match early exit:** Any local file whose full basename (including the
version suffix, e.g. `_v1.3.txt`) already exists in staging is skipped entirely. These
files are definitively unchanged — no MD5 computation, no join entry. Both sides of the
join are filtered to exclude these basenames so they cannot appear as false `REMOVED`
entries.

**Phase 2 — version-aware join:** The remaining local files (those with a new or changed
basename) are joined against the staging snapshot using a `join_key`: the basename with
the version suffix stripped.

```
human-precovid-sed-adu_t06-muscle_transcript-rna-seq_qc-norm_log-cpm_v1.3.txt
  → join_key: human-precovid-sed-adu_t06-muscle_transcript-rna-seq_qc-norm_log-cpm
```

This means a file bumped from `v1.3` to `v1.4` is matched to its staging counterpart
rather than appearing as a new `ADDED` file alongside a spurious `REMOVED`. The joined
result carries both `old_version` (from staging) and `new_version` (from local).

Change categories produced:

| Category | Definition |
|---|---|
| `ADDED` | Local file has no matching `join_key` in staging |
| `MODIFIED` | `join_key` matches in both; version bumped or content differs |

Staging files with no local counterpart (intentional removals) are excluded from the
diff output. Handle deliberate deletions manually before running Step 4.

**Output — `diffs/diff_{timestamp}.tsv`** (also copied to `diffs/latest.tsv`):

| Column | Description |
|---|---|
| `local_path` | Absolute path to the local source file |
| `gcs_path` | Derived destination GCS path for upload |
| `production_gcs_path` | GCS path of the matched file in the production bucket |
| `join_key` | Version-stripped basename used for matching |
| `ome` | Ome identifier parsed from the file name |
| `tissue_code` | BIC tissue code parsed from the file name |
| `data_category` | `qc-norm`, `metadata`, or `da` |
| `data_details` | Sub-category string (e.g. `log-cpm`, `dream-acute`) |
| `old_version` | Version string of the matched staging file |
| `new_version` | Version string of the local file |
| `old_md5` | MD5 of the staging file (from Step 1 snapshot) |
| `new_md5` | MD5 of the local file |
| `change_type` | `ADDED` or `MODIFIED` |

---

### Step 4 — Upload changed files to staging (`04_upload_to_staging.R`)

**Goal:** Apply the diff from Step 3 to the staging bucket: upload new and modified
files, and remove the old versioned file for any `MODIFIED` entry.

**What it does:**

- Reads `diffs/latest.tsv`
- For `ADDED` and `MODIFIED` entries: `gsutil cp {local_path} {gcs_path}`
- After each upload, re-runs `gsutil stat` and decodes the base64 MD5 to hex for
  comparison against the local MD5 — fails fast on any mismatch
- For `MODIFIED` entries (after verification): removes the old versioned file from
  staging by deriving its staging path from `production_gcs_path`
- Records upload success/failure, both MD5s, and old-version removal status per file
  in `logs/upload_log_{timestamp}.tsv`
- Stops if any upload fails MD5 verification or any old version cannot be removed

**Why this matters:** This is the only step that mutates the staging bucket. Verifying MD5
post-upload catches silent corruption or partial transfers before structure checks run.
Removing the old versioned file prevents two versions of the same file from co-existing
in staging.

---

### Step 5 — Validate staging folder structure (`05_validate_structure.R`)

**Goal:** Assert that the staging bucket conforms to the required directory hierarchy,
that every expected file (by ome × tissue × data_category combination) is present with
the correct column schema, and that the staged data values are consistent with the
currently installed package.

**What it does:**

**Structure checks:**

- Sources `required_structure.R` (declarative manifest: each row specifies `ome`,
  `tissue`, `tissue_code`, `data_category`, `data_details`, `gcs_subdir`, `required`)
  and `expected_columns.R` (per-ome column schemas)
- Queries `gsutil ls -R {STAGING_BUCKET}` and parses results
- For each manifest row, selects the highest-versioned matching file in staging
- Downloads each file once; checks that required columns are present and forbidden
  columns are absent
- Files matching `data_details = "removed-samples"` are treated as optional and pass
  automatically if absent
- Warns (does not fail) for subdirectories missing an HTML QC report
- Writes `logs/structure_validation_{timestamp}.tsv` — one row per manifest entry with
  `status` = `PASS | FAIL` and a `reason` column
- Stops if any row marked `required = TRUE` has status `FAIL`

**Value checks (via `qc_norm_visualization_helpers.R`):**

- Downloaded files are cached in memory (keyed by `ome__tissue__type`) during the
  structure pass — no second download needed
- For each ome × tissue with a cached `qc_norm` and `sample_metadata`, calls
  `compare_qc_norm()` which compares sample overlap, feature counts, and expression
  distributions against the currently installed `MotrpacHumanPreSuspensionData`
- Value-check messages are appended to the validation log

---

## File Structure

```
google_cloud_bucket_checks/
├── README.md                        # this document
├── config.R                         # bucket paths, version strings, local dirs
├── required_structure.R             # declarative manifest of expected GCS files
├── expected_columns.R               # per-ome expected column schemas
├── qc_norm_visualization_helpers.R  # compare_qc_norm() and read_tsv() helpers (used by Step 5)
├── run_validation_pipeline.R        # top-level driver; sources steps 1–5 in order
├── 01_snapshot_production.R
├── 02_copy_to_staging.R
├── 03_diff_local_vs_staging.R
├── 04_upload_to_staging.R
├── 05_validate_structure.R
├── snapshots/                       # output: production snapshots (auto-created)
├── diffs/                           # output: per-run diff tables (auto-created)
└── logs/                            # output: upload, validation, copy logs (auto-created)
```

---

## Running the Pipeline

### Full automated run

```r
source("data-raw/google_cloud_bucket_checks/config.R")
source("data-raw/google_cloud_bucket_checks/run_validation_pipeline.R")
```

### Step-by-step (recommended for first run or debugging)

```r
source("config.R")
source("01_snapshot_production.R")   # produces snapshots/latest.tsv
source("02_copy_to_staging.R")       # populates STAGING_BUCKET
source("03_diff_local_vs_staging.R") # produces diffs/latest.tsv
source("04_upload_to_staging.R")     # mutates STAGING_BUCKET
source("05_validate_structure.R")    # structure + value checks; produces logs/
```

### After the pipeline passes

1. Review `logs/structure_validation_*.tsv` — confirm all required rows are `PASS`
2. Review the value-check messages appended to the validation log
3. Rebuild affected `data/*.rda` objects using scripts in `data-raw/`
4. Open PRs in both downstream package repos with the updated data objects
5. After both PRs merge, promote the staging bucket to production:
   ```bash
   gsutil -m rsync -r {STAGING_BUCKET}/ {NEW_PRODUCTION_BUCKET}/
   ```

---

## Package Ownership Reference

| GCS data category | Staging subdirectory | Affected package |
|---|---|---|
| `qc-norm` | `{ome}/qc-norm/` | Data package |
| `metadata` (samples / features) | `{ome}/metadata/` | Data package |
| `DA` (differential analysis) | `{ome}/DA/` | Analysis package |
| `resources` (GMT sets, feature maps) | `resources/` | Both |
| QC reports (HTML) | `{ome}/` | Data package |

---

## Naming Convention Reference

All files follow the BIC standard enforced by `write_with_path_name()`:

```
human-precovid-sed-adu_{tissue_code}_{ome}_{data_category}_{data_details}_v{version}.txt
```

Example:
```
human-precovid-sed-adu_t06-muscle_transcript-rna-seq_qc-norm_log-cpm_v1.4.txt
```

Tissue codes are looked up from `OME_TISSUE_CODE` via `.match_ome_tissue_code()`.
