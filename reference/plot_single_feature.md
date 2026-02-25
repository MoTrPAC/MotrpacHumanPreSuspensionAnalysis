# Single feature plot function

Creates a single plot or faceted plot for a given feature or gene across
the selected tissues. Single plots are optimized for H = 1.54 in and W =
1.225 inches. If a right hand legend is added, the width should be
1.715. This can be toggled as needed by saving output and editing as a
ggplot object

## Usage

``` r
plot_single_feature(
  feature,
  selected_tissues = "all",
  selected_omes = "all",
  p_level = 0.05,
  output_file = NULL,
  scale_factor = 1,
  color_time_labels = FALSE,
  include_legend = TRUE,
  legend_position = "right",
  verbose = TRUE,
  epigen = FALSE
)
```

## Arguments

- feature:

  A single character string corresponding to a feature_id, gene_symbol,
  or refmet_name within the human feature-to-gene table

- selected_tissues:

  character; one of tissue_available_list.

- selected_omes:

  character; one of ome_available_list.

- p_level:

  Numeric threshold for adjusted p value significance in differential
  analysis (so this would make individual points highlighted in black if
  below this threshold)

- output_file:

  a file path if desired, to autoatically save output in the standard
  size of height 1.54, width 1.225, recommended if output is known to be
  a single plot, not faceted

- scale_factor:

  simple way to scale plot if larger sizes are needed for posters or
  talks. recommend integers only

- color_time_labels:

  toggle TRUE/FALSE allows color tiles corresponding to time point
  colors to be used instead of x-axis

- include_legend:

  toggle TRUE/FALSE to include a legend

- legend_position:

  if include_legend == TRUE, position can be selected

- verbose:

  logical; toggle to include verbose information

- epigen:

  logical; toggle to include epigenetic features (only atac offered for
  this function). If you include a gene name, this could result in many
  many epigenetic features mapping to the one object.

## Value

a ggplot

## Author

Daniel Katz and Christopher Jin

## Examples

``` r
if (FALSE) { # \dontrun{
plot_single_feature_test(feature="CAR 10:0",
  p_level = 0.05,
  selected_tissues = "all",
  selected_omes = "all",
  output_file = paste0("~/single_plot_test_scale.pdf"),
  scale_factor = 1,
  color_time_labels = F) +
theme(legend.position = "bottom")


plot_single_feature_test(feature="TFEB",
 p_level = 0.05,
 selected_tissues = "muscle",
 selected_omes = "prot-pr",
 output_file = paste0("~/single_plot_test_scale_tfeb_prot.png"),
 scale_factor = 1,
 color_time_labels = T,
 include_legend = F)
 } # }
```
