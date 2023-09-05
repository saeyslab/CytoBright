#' Estimate Stain Index for one detector
#'
#' @param ff        FlowFrame. Assumed to be compensated and pregated.
#' @param detector  Detector to evaluate. Should be a column name of ff.
#' @param cutoff    Threshold between positive and negative population.
#' @param meta      Dataframe of 1 row with relevant meta information.
#'                  The results will be appended to this.
#'
#' @importFrom flowCore exprs keyword parameters nrow
#' @importFrom stats quantile
#' @importFrom Biobase pData
#' @export
estimate_SI <- function(ff,
                        detector,
                        cutoff,
                        meta = data.frame(matrix(NA, nrow = 1, ncol = 0))) {
  detector <- unname(detector)

  pos <- flowCore::exprs(ff)[, detector] >= cutoff
  neg <- flowCore::exprs(ff)[, detector] < cutoff

  meta["Detector"] <- detector
  meta["Voltage"] <- as.numeric(
    flowCore::keyword(ff, paste0(
      "$P",
      which(colnames(ff) == detector),
      "V"
    ))
  )
  meta["Cutoff"] <- cutoff
  meta["Pos_count"] <- sum(pos)
  meta["Neg_count"] <- sum(neg)
  meta["MFI_pos"] <- stats::quantile(flowCore::exprs(ff)[pos, detector], 0.50)
  meta["MFI_neg"] <- stats::quantile(flowCore::exprs(ff)[neg, detector], 0.50)
  meta["Max_neg"] <- stats::quantile(flowCore::exprs(ff)[neg, detector], 0.95)
  meta["Min_neg"] <- stats::quantile(flowCore::exprs(ff)[neg, detector], 0.05)
  meta["rSD"] <- (meta["Max_neg"] - meta["Min_neg"]) / 3.29
  meta["SI"] <- (meta["MFI_pos"] - meta["MFI_neg"]) / (2 * meta["rSD"])

  pData <- Biobase::pData(flowCore::parameters(ff))
  limit <- pData[pData$name == detector, "maxRange"]
  meta["Pctg_OutOfRange"] <- sum(flowCore::exprs(ff)[, detector] >= limit) /
    flowCore::nrow(ff)

  return(meta)
}
