#' MICE imputation
#'
#' @author Natalie M Clark, Stephanie Vartany
#' @param data matrix for imputation
#' @param na_max maximum proportion of missing values allowed (default 0.4)
#' @param num_imps number of multiple imputations (default 15)
#' @param seed set random seed for consistent imputation results (default 2023)
#' @param num_cores number of cores for parallelization (default 1, no parallelization)
#'
#' @importFrom dplyr group_by
# #' @importFrom mice futuremice mice
#commenting this out because we don't want mice to be a formal package requirement
#because basically only 1-2 analysts max will actually want to run this function
#'
#' @return the imputed data matrix
run_mice <- function(data,
                     na_max=0.4,
                     num_imps=15,
                     seed=2023,
                     num_cores=1){
  # check_package_installation(pkg = "mice")

  # typically MICE wants the matrix as samples x features (feature-wise imputation), but this takes much too long.
  #We have obtained good results with sample-wise imputation which we perform here. With parallelization, feature-wise may be possible, but is likely not worth it.

  #filter using na_max
  #MICE seems to have optimal results when na_max=0.4, after that it starts to drop off
  #this is likely dataset dependent
  data_filt <- data[rowSums(is.na(data))/dim(data)[2] <= na_max,]

  # run mice, this creates a mice-specific output object
  # m is number of iterations, set to 15 (default was 5)
  # in mice v3.15 and later you can run in parallel using futuremice()
  # set seed so result is the same for the same dataset (for futuremice, use parallelseed)
  print("Imputing using MICE")
  if(num_cores>1){
    mice_out <- mice::futuremice(data_filt, m=num_imps, parallelseed=seed, n.core=num_cores, print=T)
  }else{
    mice_out <- mice::mice(data_filt, m=num_imps, seed=seed, print=T)
  }
  print("Imputation done")

  # collect the mice object (with all iterations)
  print("Aggregating imputations")
  imp_data <- mice::complete(mice_out, 'long')

  # aggregate across all iterations
  avg_imp_data <- imp_data %>%
    dplyr::group_by(.id) %>%
    summarize(across(.fns = mean)) %>%
    dplyr::select(-c('.id', '.imp'))
  avg_imp_data <- as.data.frame(avg_imp_data)
  rownames(avg_imp_data) <- rownames(data_filt)
  print("Aggregation done")

  #return the aggregated matrix
  return(avg_imp_data)
}
