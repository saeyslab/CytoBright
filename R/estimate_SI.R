#' Estimate Stain Index for one detector
#'
#' @param ff        FlowFrame. Assumed to be compensated and pregated.
#' @param detector  Detector to evaluate. Should be a column name of ff.
#' @param cutoff    Threshold between positive and negative population.
#' @param meta      Dataframe of 1 row with relevant meta information.
#'                  The results will be appended to this.
#'
#' @importFrom flowCore exprs keyword parameters
#' @importFrom stats quantile qnorm pnorm
#' @importFrom Biobase pData
#' @export
estimate_SI <- function(ff,
                        detector,
                        cutoff,
                        meta = data.frame(matrix(NA, nrow = 1, ncol = 0))) {
  detector <- unname(detector)

  exprs_ff_detector <- flowCore::exprs(ff)[, detector]
  pos <- exprs_ff_detector >= cutoff
  neg <- exprs_ff_detector < cutoff

  meta["Detector"] <- detector
  # some FCS files don't have voltage keyword
  voltage <-
    flowCore::keyword(ff, paste0(
      "$P",
      which(colnames(ff) == detector),
      "V"
    ))[[1]]
  meta["Voltage"] <- ifelse(is.null(voltage), 0, as.numeric(voltage))
  meta["Cutoff"] <- cutoff
  meta["Pos_count"] <- sum(pos)
  meta["Neg_count"] <- sum(neg)
  meta["MFI_pos"] <- stats::quantile(exprs_ff_detector[pos], 0.50)
  meta["MFI_neg"] <- stats::quantile(exprs_ff_detector[neg], 0.50)
  meta["Max_neg"] <- stats::quantile(exprs_ff_detector[neg], 0.95)
  meta["Min_neg"] <- stats::quantile(exprs_ff_detector[neg], 0.05)
  meta["rSD"] <- (meta["Max_neg"] - meta["Min_neg"]) /
    (stats::qnorm(0.95) - stats::qnorm(0.05))
  meta["SI"] <- (meta["MFI_pos"] - meta["MFI_neg"]) / (2 * meta["rSD"])

  # Out of range events
  pData <- Biobase::pData(flowCore::parameters(ff))
  limit <- pData[pData$name == detector, "maxRange"]
  meta["Pctg_OutOfRange"] <- sum(exprs_ff_detector >= limit) /
    length(exprs_ff_detector)

  return(meta)
}
