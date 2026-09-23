# PTM-SEA results

Combined PTM-SEA output (`ptm-sea-results-combined.gct`, GCT 1.3) for the prot-ph
EE-CON and RE-CON contrasts, one file per tissue. Files are unmodified copies of
the sources below; only the names changed. `PTMSEA_RESULTS.R` builds the
`PTMSEA_RESULTS` data object from them.

| File | Signatures x contrasts | Source |
|---|---|---|
| `muscle_prot-ph_ptmsea_EE-RE-vs-CON_combined.gct` | 506 x 6 | precovid-analyses `figures/muscle/Figure5/ptm-sea-results-combined.gct`, commit `fe9936e` (Natalie Clark, 2026-02-18); unchanged in PR #104 |
| `adipose_prot-ph_ptmsea_EE-RE-vs-CON_combined.gct` | 437 x 2 | Supplied by Cheehoon Ahn; identical to precovid-analyses `figures/adipose/Files/n3_ptm-sea-results-combined.gct`, commit `020da87` (2025-11-09), read by `figures/adipose/precawg_adi_pathway_v2.R` |

Git blob hashes: muscle `0a1e149bea9a30b0e8d943876f08d0903b5bec5b`, adipose
`37bf24d6d0f2168ea7e3c7fcccf8b7fb1c96d722`.

The adipose run in precovid-analyses PR #104 (`figures/adipose/Files/PTMSEA_v2.0/`)
is not used.

Neither file came with a parameters file. The smallest signature overlap in
every contrast of both files is 5 sites, consistent with `min.overlap = 5`. The
input z-matrix for each run is not recorded.
