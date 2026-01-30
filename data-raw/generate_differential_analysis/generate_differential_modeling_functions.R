#' Fit linear mixed models using the dream framework
#'
#' This function fits linear mixed models for differential analysis using the
#' \code{dream} framework from the \code{variancePartition} package. It supports
#' both acute and training analyses internally, with contrast definitions and
#' model formulas determined by the supplied metadata and analysis type.
#'
#' For acute analyses, contrasts are generated using a restricted, baseline-anchored
#' contrast set to improve statistical stability and reduce computational burden.
#' Training contrasts are generated separately when requested.
#'
#' While the function is capable of fitting both acute and training models,
#' **only acute exercise results are publicly released** due to limited sample
#' sizes and power constraints in training cohorts.
#'
#' @param expression_object
#' A model-ready expression object (e.g., \code{edgeR::DGEList} or matrix),
#' with columns corresponding to samples.
#'
#' @param acute_vs_training
#' Character scalar specifying whether the analysis corresponds to
#' \code{"acute"} or \code{"training"} exercise.
#'
#' @param process_metadata
#' A list returned by \code{process_covariates}, containing processed metadata,
#' model formulas, and design information.
#'
#' @param voom
#' Logical indicating whether voom-based precision weights should be applied
#' prior to model fitting. Typically \code{TRUE} for RNA-seq and ATAC-seq data.
#'
#' @param parallel
#' Logical indicating whether model fitting should be parallelized using
#' \code{BiocParallel}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reorders metadata to match expression columns
#'   \item Generates contrasts appropriate to the analysis type
#'   \item Constructs contrast matrices using \code{makeContrastsDream}
#'   \item Optionally applies voom-based weighting
#'   \item Fits linear mixed models using \code{dream}
#'   \item Applies empirical Bayes moderation via \code{eBayes}
#' }
#'
#' @return
#' A fitted \code{dream} model object with moderated statistics.
#'
#' @note
#' Training models are fitted for internal analyses only and are not included
#' in public data releases.
#'
#' @keywords internal
#' @author christopher jin

run_dream = function(expression_object = NULL,
                     acute_vs_training = NULL,
                     process_metadata,
                     voom = FALSE,
                     parallel = FALSE){
  # check_package_installation("variancePartition")
  if(parallel) {
    num_cores = parallel::detectCores() - 2
    param <- SnowParam(num_cores, "SOCK", progressbar = TRUE)
  }
  meta_matrix = process_metadata$metadata
  meta_matrix = meta_matrix[match(colnames(expression_object), rownames(meta_matrix)), ] #reorder so they're in the same order as colnames of raw_counts
  #in theory, ^ they should already be the same, because they're matched beforehand...but I'm keeping it just in case.
  if(acute_vs_training == "acute") contrast_expressions = .generate_contrasts_acute(meta_matrix) #generate contrasts in a consistent way
  if(acute_vs_training == "training") contrast_expressions = .generate_contrasts_training(meta_matrix)

  if(acute_vs_training == "acute")  formula = process_metadata$full_formula %>% as.formula()
  if(acute_vs_training == "training") formula = process_metadata$training_formula %>% as.formula()

  L = variancePartition::makeContrastsDream(formula,
                                            meta_matrix,
                                            contrasts = contrast_expressions)
  if (voom){
    if (parallel){
      suppressWarnings({expression_object <- variancePartition::voomWithDreamWeights(expression_object, formula, meta_matrix, BPPARAM = param)}) #only for atac, rna-seq
    }else{
      expression_object <- variancePartition::voomWithDreamWeights(expression_object, formula, meta_matrix) #only for atac, rna-seq
    }
  }

  if (parallel){
    suppressWarnings({fit = variancePartition::dream(expression_object, formula, meta_matrix, L = L, BPPARAM = param)})
  }else{
    fit = variancePartition::dream(expression_object, formula, meta_matrix, L = L)
  }
  fit = variancePartition::eBayes(fit)
  return(fit)
}


#' Convert dream model fits into a standardized differential analysis table
#'
#' This internal helper function extracts contrast-specific results from a fitted
#' \code{dream} model and converts them into a standardized, long-format table
#' suitable for downstream analysis, visualization, and public release.
#'
#' Only coefficients corresponding to explicit contrasts (i.e., containing
#' subtraction operators) are retained. For each contrast, full summary statistics
#' are extracted and annotated with model, tissue, and omics metadata.
#'
#' @param fit
#' A fitted \code{dream} model object returned by \code{run_dream}.
#'
#' @param metadata
#' Original sample metadata used for model fitting.
#'
#' @param formula
#' The model formula used to fit the mixed model.
#'
#' @param tissue
#' Character scalar specifying the tissue analyzed.
#'
#' @param ome
#' Character scalar specifying the omics layer analyzed.
#'
#' @details
#' For each retained contrast, the function extracts:
#' \itemize{
#'   \item Effect sizes and standard errors
#'   \item Moderated test statistics
#'   \item Confidence intervals
#'   \item Raw and adjusted p-values
#'   \item Model degrees of freedom and log-likelihood
#' }
#'
#' Results are concatenated across contrasts and sorted by adjusted p-value.
#'
#' @return
#' A data frame containing standardized differential analysis results across
#' all contrasts for the given tissue–omics combination.
#'
#' @note
#' This function does not perform filtering beyond contrast selection; all
#' downstream significance thresholds are applied elsewhere.
#'
#' @keywords internal output
#' @author christopher jin

.convert_dream_output = function(fit,
                                 metadata = NULL,
                                 formula = NULL,
                                 tissue = NULL,
                                 ome = NULL){
  full_contrasts = colnames(coef(fit))
  comparison_subset = full_contrasts[grep("-", full_contrasts)]
  #subset to only the models with some type of actual contrast
  res_tissue = data.frame() #final output file
  for(contrast in comparison_subset){
    res = variancePartition::topTable(fit, coef = contrast, number = Inf, p.value = 1, confint = TRUE)
    res_single = res %>%
      dplyr::mutate(degrees_of_freedom = fit$rdf) %>%
      dplyr::mutate(logLik = fit$logLik) %>%
      dplyr::mutate(feature_id = rownames(.)) %>%
      dplyr::mutate(full_model = contrast) %>%
      dplyr::mutate(assay = ome) %>%
      dplyr::mutate(contrast = contrast) %>%
      dplyr::mutate(full_model = formula) %>%
      dplyr::mutate(tissue = tissue) %>%
      dplyr::rename(p_value = P.Value) %>%
      dplyr::rename(adj_p_value = adj.P.Val) %>%
      dplyr::select(assay,
                    feature_id,
                    z.std,
                    logFC,
                    CI.L, CI.R,
                    degrees_of_freedom,
                    logLik,
                    t,
                    AveExpr,
                    p_value, adj_p_value,
                    contrast,
                    full_model)

    res_tissue = rbind(res_tissue, res_single)
  }
  res_tissue = res_tissue %>% dplyr::arrange(adj_p_value)
  return(res_tissue)
}

#' Construct covariate metadata and model formulas for mixed model analysis
#'
#' This function defines, scales, and encodes covariates for use in linear mixed
#' model differential analysis. Covariate inclusion is controlled by omics layer
#' and tissue, with optional exclusion of technical covariates to support
#' sensitivity analyses.
#'
#' The function returns a structured list containing processed metadata matrices,
#' model formulas for acute and training analyses, and documentation of included
#' technical and design covariates.
#'
#' @param meta
#' A data frame of sample metadata containing participant identifiers, group
#' assignments, timepoints, and candidate covariates.
#'
#' @param selected_ome
#' Character scalar specifying the omics layer being analyzed.
#'
#' @param tissue_input
#' Character scalar specifying the tissue being analyzed.
#'
#' @param include_technical
#' Logical indicating whether technical covariates should be included in the
#' mixed model. Default is \code{TRUE}.
#'
#' @param custom_covariates
#' Optional data frame specifying custom covariate definitions. If \code{NULL},
#' a package-defined covariate configuration file is used. See methods for hose these
#' covariate configurations were originally derived.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Selects covariates relevant to the chosen ome and tissue
#'   \item Scales numerical covariates
#'   \item Encodes categorical covariates as factors
#'   \item Constructs interaction terms for group and timepoint
#'   \item Assembles fixed and random effect model formulas
#' }
#'
#' Acute analyses use participant-level random intercepts, while training
#' analyses include visit-level random slopes.
#'
#' @return
#' A named list containing:
#' \itemize{
#'   \item Processed metadata matrix
#'   \item Acute and training model formulas
#'   \item Lists of technical and design covariates
#'   \item Non–mixed model formula for diagnostic use
#' }
#'
#' @note
#' Although both acute and training formulas are generated, only acute analyses
#' are intended for public release due to training sample size limitations.
#'
#' @keywords modeling covariates
#' @author christopher jin

process_covariates = function(meta,
                              selected_ome,
                              tissue_input,
                              include_technical = TRUE,
                              custom_covariates = NULL){
  covariates_return = list() #output list for the end part
  covariates_return[["original_meta"]] = meta

  if(!is.null(custom_covariates)){
    input_covariates = custom_covariates
  }else{
    input_covariates = COVARIATES_FILE
  }
  covariates = input_covariates %>%
    as.data.frame() %>%
    dplyr::filter(ome == selected_ome) %>%
    dplyr::filter(tissue == 'all' | tissue == tissue_input)

  num_cov = covariates %>% dplyr::filter(data_type == "numerical") #numerical covariates
  factor_cov = covariates %>% dplyr::filter(data_type == "factor")

  sel_meta = meta %>%
    dplyr::select(all_of(covariates$covariate)) %>%
    dplyr::mutate(across(all_of(num_cov$covariate), ~ scale(.) %>% as.numeric())) %>%
    dplyr::mutate(across(all_of(factor_cov$covariate), ~ as.factor(.) %>% droplevels())) %>%
    dplyr::mutate(group_timepoint = droplevels(interaction(randomGroupCode, Timepoint))) %>%
    dplyr::mutate(visit_group_timepoint = droplevels(interaction(visitcode, randomGroupCode, Timepoint)))

  technical_covs = covariates %>% filter(tech_or_design == "Technical")
  full_formula = names(sel_meta)[!names(sel_meta) %in% c("randomGroupCode", "Timepoint", "visitcode", "pid", "group_timepoint", "visit_group_timepoint")] #remove these from the character vector
  #the purpose of the design covariates section is to make a model.matrix()
  design_covs = c(full_formula[!full_formula %in% as.character(technical_covs$covariate)], "group_timepoint")
  if(!include_technical){ #remove for any modeling where some covariates have been regressed out
    full_formula = full_formula[!full_formula %in% technical_covs$covariate]
  }
  formula_string = paste(full_formula, collapse = " + ")
  #we readd group_timepoint first because of the way some contrast matrixes drop values in case of interactions w other levels of factors in the contrast matrixes
  formula_string_full = paste("~ 0 + group_timepoint + ", formula_string, "+ (1 | pid)")
  formula_string_training = paste("~ 0 + visit_group_timepoint + ", formula_string, "+ (visitcode | pid)")

  non_mixed_model = paste("~ 0 + group_timepoint + ", formula_string)

  #so i remove group_timepoint above and then make sure that it comes first because the order of the string can sometimes
  #actually change the contrast matrix formed and which columns are dropped in terms of the contrast comparisons.
  covariates_return[["technical_cov"]] = technical_covs
  covariates_return[["design_cov"]] = design_covs

  covariates_return[["full_formula"]] = formula_string_full
  covariates_return[["training_formula"]] = formula_string_training

  covariates_return[["metadata"]] = sel_meta #so this is with all the tech/num cov in the correct format
  covariates_return[["non_mixed_model"]] = non_mixed_model

  return(covariates_return)
}


#' Generate a restricted set of acute exercise contrasts for linear mixed models
#'
#' This internal helper function constructs a curated set of contrast expressions
#' for acute exercise analyses, designed to capture biologically interpretable
#' within-group timepoint effects and between-group differences while avoiding
#' the combinatorial explosion associated with fully pairwise contrast generation.
#'
#' Rather than enumerating all possible timepoint × group comparisons, which
#' substantially increases computational burden and multiple-testing penalties,
#' this function defines a semi-manual contrast scheme anchored to the
#' \code{pre_exercise} baseline. This approach balances interpretability,
#' statistical efficiency, and computational feasibility.
#'
#' Only acute exercise contrasts are generated. Training-related contrasts are
#' intentionally excluded, as only acute analyses are publicly released due to
#' sample size and power considerations.
#'
#' @param metadata
#' A data frame of sample metadata containing at least a \code{Timepoint} column.
#' Timepoints must include \code{pre_exercise} and one or more post-baseline
#' acute exercise timepoints.
#'
#' @details
#' For each post-baseline acute timepoint, the following contrasts are generated
#' when applicable:
#' \itemize{
#'   \item Endurance vs Control (change from pre-exercise)
#'   \item Endurance within-group change from pre-exercise
#'   \item Resistance vs Control (change from pre-exercise)
#'   \item Resistance within-group change from pre-exercise
#'   \item Endurance vs Resistance (change from pre-exercise)
#'   \item Control within-group change from pre-exercise
#' }
#'
#' For timepoints lacking resistance blood draws (e.g.,
#' \code{during_20_min}, \code{during_40_min}), resistance-related contrasts
#' are automatically omitted.
#'
#' In addition, baseline (pre-exercise) group contrasts are included to
#' characterize between-group differences prior to exercise onset.
#'
#' @return
#' A character vector of contrast expressions suitable for use with
#' \code{limma}- or \code{dream}-based modeling frameworks.
#'
#' @note
#' This restricted contrast set is intentionally conservative and is used
#' to support stable inference in settings with limited sample sizes.
#' Fully pairwise contrast generation is deliberately avoided.
#'
#' @keywords internal

.generate_contrasts_acute = function(metadata){
  pre_contrast_expressions <- c()
  timepoints = unique(metadata$Timepoint)
  for(tp in timepoints){
    # message(tp)
    if (!tp == 'pre_exercise'){
      # meta_tp = metadata %>% filter(Timepoint == tp)
      contrast_Endur_Cntrl = sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise - group_timepointADUControl.%s + group_timepointADUControl.pre_exercise", tp, tp)
      contrast_Endur = sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise", tp)
      contrast_Resist_Cntrl = sprintf("group_timepointADUResist.%s - group_timepointADUResist.pre_exercise - group_timepointADUControl.%s + group_timepointADUControl.pre_exercise", tp, tp)
      contrast_Resist = sprintf("group_timepointADUResist.%s - group_timepointADUResist.pre_exercise", tp)
      contrast_Endur_Resist = sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise - group_timepointADUResist.%s + group_timepointADUResist.pre_exercise", tp, tp)
      contrast_Cntrls = sprintf("group_timepointADUControl.%s - group_timepointADUControl.pre_exercise", tp)
      if (tp == 'during_20_min' | tp == 'during_40_min') {contrast_Resist_Cntrl = NULL; contrast_Endur_Resist = NULL; contrast_Resist = NULL} #resistance group doesn't get blood draws here
      pre_contrast_expressions = c(pre_contrast_expressions,
                                   contrast_Endur_Cntrl,
                                   contrast_Endur,
                                   contrast_Resist_Cntrl,
                                   contrast_Resist,
                                   contrast_Endur_Resist,
                                   contrast_Cntrls)
    }
  }

  pre_ex_endur_res = "group_timepointADUEndur.pre_exercise - group_timepointADUResist.pre_exercise"
  pre_ex_res_cntrl = "group_timepointADUResist.pre_exercise - group_timepointADUControl.pre_exercise"
  pre_ex_endur_cntrl = "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise"
  pre_contrast_expressions = c(pre_contrast_expressions, pre_ex_endur_res, pre_ex_res_cntrl, pre_ex_endur_cntrl)
  return(pre_contrast_expressions)
}


