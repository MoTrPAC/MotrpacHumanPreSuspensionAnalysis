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
  repo_local_dir = NULL,
  verbose = TRUE,
  epigen = FALSE,
  qc_data = NULL
)
```

## Arguments

- feature:

  A single character string corresponding to a feature_id, gene_symbol,
  or refmet_name within the human feature-to-gene table

- selected_tissues:

  character; one of tissue_available_list.

- selected_omes:

  character; one of ome_available_list. The clinical chemistry omes in
  [`clinical_ome_list()`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/clinical_ome_list.md)
  are plotted only when they are requested, either by name or through
  `"all"`; asking for another ome no longer returns clinical chemistry
  alongside it. `"metab"` does not imply `"metab-t-clinical"`.

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

- repo_local_dir:

  Deprecated and ignored; passed through to
  [`load_differential_analysis`](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/reference/load_differential_analysis.md),
  which prints a message if it is supplied.

- verbose:

  logical; toggle to include verbose information

- epigen:

  logical; toggle to include epigenetic features (only atac offered for
  this function). If you include a gene name, this could result in many
  many epigenetic features mapping to the one object. Only the
  epigenomic files for the requested tissues and omes are downloaded.

- qc_data:

  optional; the nested list returned by
  `MotrpacHumanPreSuspensionData::load_qc()`, for users with access to
  the individual-level data. The shipped `*_SUM_STATS` objects cover
  only the features in the differential analysis, and for the epigenomic
  omes only the features with `adj_p_value < 0.05` in at least one
  contrast. When a requested feature has differential analysis but no
  summary statistics, they are computed from `qc_data` the way the
  shipped objects are: `qc_norm` restricted to `visitcode == "ADU_BAS"`
  samples, summarised per `randomGroupCode` and `Timepoint`. Load it
  with the tissues and omes being plotted, and `epigen = TRUE` for the
  epigenomic omes.

## Value

a ggplot

## Features without summary statistics

The lines, points and error bars are drawn from the summary statistics,
not from the differential analysis. A feature that is in the
differential analysis but not in the summary statistics, most often an
epigenomic feature that is not significant in any contrast, has nothing
to draw. The function then says so with a
[`message()`](https://rdrr.io/r/base/message.html) naming the features
and returns a plot with empty panels for them, unless `qc_data` is
supplied to compute the missing statistics.

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
