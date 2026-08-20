#' Identify channel with highest value at a certain quantile level
#'
#' @param ff FlowFrame
#' @param q quantile for which the values should be evaluated
#'
#' @importFrom flowCore exprs colnames
#' @importFrom stats quantile
#' @export
get_highest_channel <- function(ff, q = 0.9, ref_ff = NULL){
  area_columns <- grep("-A", colnames(flowCore::exprs(ff)), value = TRUE)
  fluor_columns <- grep("SC|LightLoss|Img", area_columns,
                        value = TRUE, invert = TRUE)

  if(!is.null(ref_ff)){
    ref <- apply(flowCore::exprs(ref_ff)[,fluor_columns], 2, quantile, q)
  } else {
    ref <- rep(0, length(fluor_columns))
  }

  signal <- apply(flowCore::exprs(ff)[,fluor_columns], 2, quantile, q)


  names(which.max(signal - ref))
}
