#' @title Pre-processing for PTM-SEA
#'
#' @description This code produces a GCT file of the selected muscle contrasts
#' which is used as input to the PTM-SEA docker.
#'
#' @param selected_tissues character; tissue to use, either adipose or muscle.
#' Default is muscle. Only one tissue may be processed at a time.
#' @param selected_contrast_type character; one or more contrast_types to use. Default is
#' \code{"exercise_with_controls"}.
#' @param selected_contrast_category character; one or more contrast_categories to use.
#' Default is \code{"c("EE-CON","RE-CON")"}
#'
#' @returns A GCT object with the properly formatted DA results. This GCT object
#' can be saved locally for use with PTM-SEA.
#'
#' @author Natalie M Clark
#'
#' @importFrom dplyr %>% select rename full_join bind_rows mutate across
#'   left_join arrange where desc everything
#'
#' @export preprocess_PTMSEA
#'
#' @examples
#' \dontrun{
#' # Produce GCT file for muscle phospho
#' gct <- preprocess_PTMSEA(selected_tissues="muscle")
#'
#' # If you want to save the file locally, use write_gct
#' write_gct(gct,"muscle_PTMSEA.gct",appenddim=F)
#' }

preprocess_PTMSEA <- function(selected_tissues = "muscle",
                              selected_contrast_type = "exercise_with_controls",
                              selected_contrast_category = c("EE-CON","RE-CON"))
{
  # check cmapR is installed
  check_package_installation(pkg = "cmapR", fun = "GCT")

  #obtain the correct DA results and filter accordingly
  if(selected_tissues=="muscle"){
    prot_da <- MotrpacHumanPreSuspensionAnalysis::MUSCLE_PROT_PH_DA %>%
      dplyr::filter(contrast_type %in% selected_contrast_type &
                      contrast_category %in% selected_contrast_category)
  }else{
    prot_da <- MotrpacHumanPreSuspensionAnalysis::ADIPOSE_PROT_PH_DA %>%
      dplyr::filter(contrast_type %in% selected_contrast_type &
                      contrast_category %in% selected_contrast_category)
  }
  #pivot wider
  prot_da_wide <- pivot_wider(prot_da,
                              id_cols="feature_id",
                              names_from="contrast_short",
                              values_from=c("logFC","z.std","p_value","adj_p_value"))

  #add feature metadata
  if(selected_tissues=="muscle"){
    prot_meta <- MUSCLE_PROT_PH_QC$feature_metadata %>% dplyr::rename(feature_id=id)
  }else{
    prot_meta <- ADIPOSE_PROT_PH_QC$feature_metadata %>% dplyr::rename(feature_id=id)
  }
  prot_da_wide <- right_join(prot_meta,prot_da_wide,by="feature_id")

  #create GCT object for PTMSEA
  #filter for fully localized sites
  prot_da_wide <- prot_da_wide[prot_da_wide$confident_site,]

  #use z-scores for the values
  rdesc <- prot_da_wide %>% dplyr::select(!starts_with('z.std'))
  mat <- prot_da_wide %>% dplyr::select(starts_with('z.std'))
  rownames(mat) <- rdesc$feature_id

  #add tissue to the contrast names
  colnames(mat) <- paste(selected_tissues,colnames(mat),sep=".")

  #create and return the GCT object
  gct <- GCT(mat=as.matrix(mat),
             rdesc=as.data.frame(rdesc),
             rid=rownames(mat),
             cid=colnames(mat))
  return(gct)
}
