# Fuzzy c-means clustering and enrichment analysis results

Results derived from fuzzy c-means (FCM) clustering of molecular
features, together with downstream gene set enrichment analyses. These
datasets summarize clustering assignments and cluster-level functional
enrichment using both CAMERA-based competitive gene set testing and
classical over-representation analysis (ORA).

The FCM results provide soft cluster membership estimates, allowing
features to partially belong to multiple clusters, while the enrichment
results characterize the biological functions associated with each
cluster.

## Usage

``` r
FCM_CLUSTERS
FCM_CAMERA
FCM_ORA
```

## Format

Objects of class `data.frame`.

- `FCM_CLUSTERS`:

  Feature-level fuzzy c-means clustering results. Typically includes
  feature identifiers, cluster centers, and membership values indicating
  the degree of association between each feature and each cluster.

- `FCM_CAMERA`:

  Cluster-level functional enrichment results obtained using CAMERA.
  Contains statistics such as enrichment direction, test statistics,
  p-values, and multiple-testing–adjusted p-values for each gene set
  evaluated within each cluster.

- `FCM_ORA`:

  Cluster-level over-representation analysis results. Summarizes gene
  sets that are significantly over-represented among features assigned
  to each cluster, along with corresponding enrichment statistics.

## Details

These datasets are intended to be used jointly. Fuzzy c-means clustering
defines coherent temporal or condition-specific feature patterns, while
CAMERA and ORA provide complementary perspectives on the functional
interpretation of those clusters.
