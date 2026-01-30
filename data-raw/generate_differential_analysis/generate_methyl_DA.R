#' Briefly, methylcap is processed separately using a different pipeline and so there is no function in this package to generate the "QC-Norm" dataset.
#' Differential analysis is performed with a GLMM to account for both the measured methylated and unmethylated counts.

.generate_methyl_inputs = function(){
  message("Generating a methylcap DA requires a seperate pipeline, see QC Notebook on the precovid-analyses repo for more info")
}
