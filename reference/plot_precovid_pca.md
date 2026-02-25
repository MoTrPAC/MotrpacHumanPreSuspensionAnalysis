# Principal component analysis for the pre-cawg group

Basically just a wrapper function to plot with shapes and colors that
fits the pre-cawg figure guidelines#'

## Usage

``` r
plot_precovid_pca(
  data_to_pca,
  sample_metadata,
  to_scale = TRUE,
  custom_title = ""
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

## Value

a ggplot object with the desired PCA plotting aesthetic.
