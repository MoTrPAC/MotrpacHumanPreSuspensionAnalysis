# @title SCION implementation for preCAWG. Combines preprocessing and network inference into one function

This function is based on the source code
`https://github.com/nmclark2/SCION`

## Usage

``` r
run_SCION(
  randomGroupCode = c("ADUResist", "ADUEndur"),
  regulators,
  targets,
  permute = NULL,
  dim = "col",
  cluster = TRUE,
  dir.name = "exp",
  weightthreshold = 0,
  normalize = FALSE,
  num.cores = 1,
  connect.hubs = TRUE,
  verbose = TRUE
)
```

## Arguments

- randomGroupCode:

  - regulators, targets, and clustering is all done on a per-modality
    basis. All the clustering is done using z-scores from the
    differential analysis

- regulators:

  the set of hypothesized tissue/omes that regulate the targets matrix,
  from .load_scion_matrixes

- targets:

  the set of hypothesized tissue/omes that are targetted by the
  regulators matrix, from `.load_scion_matrixes`

- permute:

  number of random permutations to perform for edge trimming. Default
  NULL, which means no permutations will be performed

- dim:

  dimension on which to permute. options are "row" or "col". default
  "col"

- cluster:

  boolean option to cluster using c-means prior to network inference.
  default TRUE

- dir.name:

  name of directory to save results. default "exp". if using
  permutations, this parameter is ignored, and the directories are named
  after the permutation number (1,2,3...)

- weightthreshold:

  threshold for edge trimming. all edges with weight \< weightthreshold
  are removed from the network. minimum value of 0. default 0.

- normalize:

  boolean to normalize edge weights to a 0-1 scale. default FALSE

- num.cores:

  number of cores to use for parallelization. num.cores-1 will be used
  for parallelization. if num.cores \< 3, no parallelization is
  performed. default 1.

- connect.hubs:

  boolean to connect the hubs between clusters. this parameter is
  ignored if clustering is not performed. default TRUE

- verbose:

  boolean to display detailed output. default FALSE

## Author

Natalie M Clark, Christopher Jin
