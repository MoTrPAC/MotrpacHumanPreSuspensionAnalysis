# canCor analysis cor correlation to principal components

Performs PCA and correlates PCs to selected metadata using
[`variancePartition::canCorPairs`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/canCorPairs.md).

## Usage

``` r
plot_precovid_cca(
  data_to_pca,
  sample_metadata,
  to_scale = TRUE,
  custom_title = "",
  interesting_variables = c("pid", "randomGroupCode", "Timepoint", "Sex", "BMI",
    "calculatedAge"),
  num_pcs = 5
)
```

## Arguments

- data_to_pca:

  The actual numerical data that should be included in analysis. must be
  able to be converted into a numerical matrix

- sample_metadata:

  The correponding metadata to `data_to_pca`. The rownames of this
  matrix must correspond to colnames of the data matrix.

- to_scale:

  logical; Whether to scale the data when performing PCA.

- custom_title:

  character; whatever title is desired for the figure

- interesting_variables:

  character; a vector of variables that are going to be correlated using
  canCors to the PCs. The \# of PCs in `num_pcs` are included by
  default.

- num_pcs:

  numeric; the number of principal components desired to be included.

## Value

a pheatmap object with the desired correlations to principal components

## Author

Christopher Jin
