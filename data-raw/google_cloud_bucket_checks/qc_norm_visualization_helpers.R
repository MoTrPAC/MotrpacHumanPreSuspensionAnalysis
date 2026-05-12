read_tsv = function(path) read.table(path, header = TRUE, sep = "\t", check.names = FALSE, quote = "")

path_for = function(local_path, ome, tissue, category, details, version = "1.2") {
  write_with_path_name(local_path = local_path, ome = ome,
                       tissue = tissue, data_category = category,
                       data_details = details, version = version,
                       return_name_only = TRUE)
}

load_regenerated_qc = function(ome, qc_norm_path, metadata_path,
                              tissues = c("muscle", "blood", "adipose"),
                              qc_norm_details = "log-cpm",
                              version = "1.2") {
  read_with_context = function(path, label) {
    tryCatch(
      read_tsv(path),
      error = function(e) stop(label, " [", path, "]: ", conditionMessage(e), call. = FALSE)
    )
  }

  result = list()
  for(tissue in tissues) {
    feat_path = path_for(metadata_path, ome, tissue, "metadata", "features", version)
    samp_path = path_for(metadata_path, ome, tissue, "metadata", "samples", version)
    qc_path = path_for(qc_norm_path, ome, tissue, "qc-norm", qc_norm_details, version)
    feature_metadata = read_with_context(feat_path, paste0(ome, "/", tissue, " feature_metadata"))
    sample_metadata = read_with_context(samp_path, paste0(ome, "/", tissue, " sample_metadata"))
    qc_norm = read_with_context(qc_path, paste0(ome, "/", tissue, " qc_norm"))
    result[[tissue]] = list(
      feature_metadata = feature_metadata,
      sample_metadata = sample_metadata,
      qc_norm = qc_norm
    )
  }
  result
}

#' Compare QC-normalized outputs between the released package and regenerated files
#'
#' Runs four sequential checks per tissue and prints a message-based summary:
#'
#' 1. qc_norm feature/sample overlap: compares row identifiers (pkg rownames vs
#'    regen feature_id column) and column identifiers (pkg colnames vs regen
#'    sample columns), listing up to 5 examples of any asymmetric entries.
#'
#' 2. qc_norm value correlation: Pearson r across the vectorized intersection of
#'    shared features x shared samples, skipping NAs. Reports r=1 as "perfectly
#'    correlated", otherwise reports r rounded to 4 decimal places with the
#'    shared dimensions.
#'
#' 3. feature_metadata ID overlap: compares pkg feature_metadata[[feature_key]]
#'    against regen feature_metadata$feature_id (regenerated files always use
#'    "feature_id"; pkg may differ, e.g. "protein_id" for proteomics). Lists up
#'    to 5 examples per side when differences exist. Reports "identical" once
#'    across all tissues if none.
#'
#' 4. sample_metadata vialLabel overlap and column value comparison: checks
#'    vialLabel identity across tissues, then for each tissue compares values in
#'    all shared columns for shared samples. Only reports columns where at least
#'    one value differs; silent if all shared column values are identical.
#'
#' @param ome ome string, e.g. "transcript-rna-seq" or "prot-pr"
#' @param existing_pkg nested list from MotrpacHumanPreSuspensionData::load_qc(),
#'   structured as tissue -> ome -> list(feature_metadata, sample_metadata, qc_norm)
#' @param regenerated nested list from load_regenerated_qc(), same structure
#' @param tissues character vector of tissues to compare
#' @param feature_key column name in feature_metadata used as the feature identifier;
#'   "feature_id" for transcriptomics/ATAC, "protein_id" for proteomics
compare_qc_norm = function(ome, existing_pkg, regenerated,
                          tissues = c("muscle", "blood", "adipose"),
                          feature_key = "feature_id") {
  message("=== ", ome, " ===")

  # 1. qc_norm: feature and sample set overlap
  message("-- qc_norm features and samples --")
  for (tissue in tissues) {
    pkg = existing_pkg[[tissue]][[ome]]
    regen = regenerated[[tissue]]

    pkg_features = rownames(pkg$qc_norm)
    regen_features = regen$qc_norm$feature_id
    pkg_samples = colnames(pkg$qc_norm)
    regen_samples = setdiff(colnames(regen$qc_norm), "feature_id")

    only_pkg_feat = setdiff(pkg_features, regen_features)
    only_regen_feat = setdiff(regen_features, pkg_features)
    only_pkg_samp = setdiff(pkg_samples, regen_samples)
    only_regen_samp = setdiff(regen_samples, pkg_samples)

    if (length(only_pkg_feat) == 0 && length(only_regen_feat) == 0) {
      message(tissue, " features: identical (n=", length(pkg_features), ")")
    } else {
      message(tissue, " features: ", length(only_pkg_feat), " only in pkg, ",
              length(only_regen_feat), " only in regen, ",
              length(intersect(pkg_features, regen_features)), " shared")
      if (length(only_pkg_feat) > 0)
        message("  only in pkg: ", paste(head(only_pkg_feat, 5), collapse = ", "),
                if (length(only_pkg_feat) > 5) paste0(" ... +", length(only_pkg_feat) - 5, " more"))
      if (length(only_regen_feat) > 0)
        message("  only in regen: ", paste(head(only_regen_feat, 5), collapse = ", "),
                if (length(only_regen_feat) > 5) paste0(" ... +", length(only_regen_feat) - 5, " more"))
    }

    if (length(only_pkg_samp) == 0 && length(only_regen_samp) == 0) {
      message(tissue, " samples: identical (n=", length(pkg_samples), ")")
    } else {
      message(tissue, " samples: ", length(only_pkg_samp), " only in pkg, ",
              length(only_regen_samp), " only in regen, ",
              length(intersect(pkg_samples, regen_samples)), " shared")
      if (length(only_pkg_samp) > 0)
        message("  only in pkg: ", paste(head(only_pkg_samp, 5), collapse = ", "),
                if (length(only_pkg_samp) > 5) paste0(" ... +", length(only_pkg_samp) - 5, " more"))
      if (length(only_regen_samp) > 0)
        message("  only in regen: ", paste(head(only_regen_samp, 5), collapse = ", "),
                if (length(only_regen_samp) > 5) paste0(" ... +", length(only_regen_samp) - 5, " more"))
    }
  }

  # 2. qc_norm: Pearson correlation across shared features x shared samples
  message("-- qc_norm value correlation (shared features x shared samples) --")
  for (tissue in tissues) {
    pkg = existing_pkg[[tissue]][[ome]]
    regen = regenerated[[tissue]]

    shared_features = intersect(rownames(pkg$qc_norm), regen$qc_norm$feature_id)
    shared_samples = intersect(
      colnames(pkg$qc_norm),
      setdiff(colnames(regen$qc_norm), "feature_id")
    )

    pkg_mat = as.matrix(pkg$qc_norm[shared_features, shared_samples])
    regen_mat = as.matrix(regen$qc_norm[match(shared_features, regen$qc_norm$feature_id), shared_samples])
    r = cor(as.vector(pkg_mat), as.vector(regen_mat), use = "complete.obs")

    if (isTRUE(all.equal(r, 1))) {
      message(tissue, ": perfectly correlated (r=1)")
    } else {
      message(tissue, ": r=", round(r, 4), " across ",
              length(shared_features), " features x ", length(shared_samples), " samples")
    }
  }

  # 3. feature_metadata: feature ID set overlap keyed on feature_key
  message("-- feature_metadata (", feature_key, ") --")
  all_feat_same = TRUE
  for (tissue in tissues) {
    pkg_feat = existing_pkg[[tissue]][[ome]]$feature_metadata[[feature_key]]
    regen_feat = regenerated[[tissue]]$feature_metadata$feature_id

    only_pkg = setdiff(pkg_feat, regen_feat)
    only_regen = setdiff(regen_feat, pkg_feat)

    if (length(only_pkg) > 0 || length(only_regen) > 0) {
      all_feat_same = FALSE
      message(tissue, ": ", length(only_pkg), " only in pkg, ",
              length(only_regen), " only in regen, ",
              length(intersect(pkg_feat, regen_feat)), " shared")
      if (length(only_pkg) > 0)
        message("  only in pkg: ", paste(head(only_pkg, 5), collapse = ", "),
                if (length(only_pkg) > 5) paste0(" ... +", length(only_pkg) - 5, " more"))
      if (length(only_regen) > 0)
        message("  only in regen: ", paste(head(only_regen, 5), collapse = ", "),
                if (length(only_regen) > 5) paste0(" ... +", length(only_regen) - 5, " more"))
    }
  }
  if (all_feat_same) message("All tissues: ", feature_key, " identical across pkg and regen")

  # 4. sample_metadata: vialLabel overlap, then value comparison for shared columns
  message("-- sample_metadata vialLabels --")
  all_samp_same = TRUE
  for (tissue in tissues) {
    pkg_samp = existing_pkg[[tissue]][[ome]]$sample_metadata$vialLabel
    regen_samp = regenerated[[tissue]]$sample_metadata$vialLabel

    only_pkg = setdiff(pkg_samp, regen_samp)
    only_regen = setdiff(regen_samp, pkg_samp)

    if (length(only_pkg) > 0 || length(only_regen) > 0) {
      all_samp_same = FALSE
      message(tissue, ": ", length(only_pkg), " only in pkg, ",
              length(only_regen), " only in regen, ",
              length(intersect(pkg_samp, regen_samp)), " shared")
    }
  }
  if (all_samp_same) message("All tissues: vialLabels identical across pkg and regen")

  message("-- sample_metadata column values (shared columns, shared samples) --")
  for (tissue in tissues) {
    pkg_sm = existing_pkg[[tissue]][[ome]]$sample_metadata
    regen_sm = regenerated[[tissue]]$sample_metadata

    shared_samples = intersect(pkg_sm$vialLabel, regen_sm$vialLabel)
    shared_cols = setdiff(intersect(colnames(pkg_sm), colnames(regen_sm)), "vialLabel")

    pkg_sub = pkg_sm[pkg_sm$vialLabel %in% shared_samples, shared_cols, drop = FALSE]
    regen_sub = regen_sm[regen_sm$vialLabel %in% shared_samples, shared_cols, drop = FALSE]
    sort_order_pkg = order(pkg_sm$vialLabel[pkg_sm$vialLabel %in% shared_samples])
    sort_order_regen = order(regen_sm$vialLabel[regen_sm$vialLabel %in% shared_samples])
    pkg_sub = pkg_sub[sort_order_pkg, , drop = FALSE]
    regen_sub = regen_sub[sort_order_regen, , drop = FALSE]
    rownames(pkg_sub) = rownames(regen_sub) = NULL

    differing = Filter(function(col) !identical(pkg_sub[[col]], regen_sub[[col]]), shared_cols)

    any_diff = FALSE
    for (col in differing) {
      n_diff = sum(pkg_sub[[col]] != regen_sub[[col]], na.rm = TRUE)
      if (n_diff > 0) {
        message(tissue, " [", col, "]: ", n_diff, " of ", length(shared_samples), " values differ")
        any_diff = TRUE
      }
    }
    if (!any_diff) message(tissue, ": all shared column values are identical")
  }

  invisible(NULL)
}
