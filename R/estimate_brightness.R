#' Estimate SI brightness from single stains
#'
#' @param single_stains Dataframe with at least the following columns:
#'                      ID, Fluorochrome, Detector, File.
#'                      ID is assumed to be unique.
#'                      Detector should correspond with the colnames of the
#'                      flow frames.
#'                      File should be full path to the single stain fcs files,
#'                      which have been preprocessed (e.g. clean up gates been
#'                      applied).
#' @param return_cells  Boolean. If true, a list with subsampled cell values
#'                      are returned, e.g. for plotting.
#' @param seed          Seed for reproducability
#' @param comp          Compensation matrix.
#' @param transform Transform parameter to be passed to find_cutoff_flowDensity.
#'                  Can either be a logical value or a transformList.
#'                  If FALSE, no transform is applied. If TRUE (default),
#'                  flowcore::estimateLogicle is called, and if this fails, the
#'                  default logicleTransform() is applied. If a transformList,
#'                  this transformList is applied. Note that the value is
#'                  returned in the original space.
#' @param unstained     If this parameter is provided, it should contain a path
#'                      to a preprocessed unstained fcs file. The positivity
#'                      threshold will be defined on this sample, and the cells
#'                      will be concatened to the single stain for the calculations.
#'
#' @importFrom flowCore read.FCS transform transformList arcsinhTransform colnames compensate
#' @importFrom flowDensity deGate
#' @importFrom PeacoQC RemoveMargins RemoveDoublets
#' @importFrom FlowSOM AggregateFlowFrames
#'
#' @export
estimate_brightness <- function(single_stains,
                                return_cells = TRUE,
                                seed = 1,
                                comp = NULL,
                                transform = TRUE,
                                unstained = NULL,
                                remove_neg_from_single = FALSE,
                                estimate_spillover = FALSE,
                                estimate_spread = FALSE,
                                silent = TRUE) {
  rownames(single_stains) <- single_stains$ID
  detectors <- unique(single_stains$Detector)


  pb <- utils::txtProgressBar(min = 0, max = nrow(single_stains), style = 3)

  values_of_interest <- c(
    "Cutoff", "MFI_pos", "MFI_neg",
    "q05_neg", "q95_neg", "rSD",
    "SI", "Voltage"
  )

  if (!is.null(unstained)) {
    ff_unstained <- flowCore::read.FCS(unstained,
                                       truncate_max_range = FALSE,
                                       emptyValue = FALSE)
  }

  if(estimate_spillover | estimate_spread){
    values_for_comp <- c(paste0("MFI_pos_", detectors),
                         paste0("MFI_neg_", detectors),
                         paste0("Comp_", detectors))
  } else {
    values_for_comp <- c()
  }

  if(estimate_spread){
    values_for_spread <- c(paste0("MFI_pos_c_", detectors),
                           paste0("MFI_neg_c_", detectors),
                           paste0("Spread_", detectors),
                           paste0("SSI_", detectors))
  } else {
    values_for_spread <- c()
  }

  SI <- data.frame(
    matrix(
      nrow = 0,
      ncol = ncol(single_stains) +
        length(values_of_interest) +
        length(values_for_comp) +
        length(values_for_spread),
      dimnames = list(
        NULL,
        c(
          colnames(single_stains),
          values_of_interest,
          values_for_comp,
          values_for_spread
        )
      )
    ),
    check.names = FALSE
  )

  cells <- data.frame(matrix(NA,
    nrow = 0,
    ncol = 2,
    dimnames = list(NULL, c("ID", "Value"))
  ))

  for(i in seq_len(nrow(single_stains))) {

    utils::setTxtProgressBar(pb, i)

    id <- single_stains$ID[i]
    file <- single_stains$File[i]
    fluor <- single_stains[i, "Fluorochrome"]
    detector <- single_stains[i, "Detector"]

    ff <- flowCore::read.FCS(file,
                             truncate_max_range = FALSE,
                             emptyValue = FALSE)

    SI[id, colnames(single_stains)] <- single_stains[id, ]


    if (!is.null(unstained)) {
      cutoff <- find_cutoff_FMO(ff_unstained, detector, 0.995)

      if(remove_neg_from_single){
        cutoff_tmp <- find_cutoff_flowDensity(ff,
                                              detector = detector,
                                              transform = transform)
        flowCore::exprs(ff) <- rbind(
          flowCore::exprs(ff)[exprs(ff)[,detector] > cutoff_tmp, ],
          flowCore::exprs(ff_unstained)) # To be further optimized

      } else {
        flowCore::exprs(ff) <- rbind(
          flowCore::exprs(ff),
          flowCore::exprs(ff_unstained)) # To be further optimized
      }

    } else {
      cutoff <- find_cutoff_flowDensity(ff,
                                        detector = detector,
                                        transform = transform)
    }

    SI_tmp <- estimate_SI(
      ff = ff,
      detector = detector,
      cutoff = cutoff
    )

    SI[id, colnames(SI_tmp)] <- SI_tmp


    if(estimate_spillover){
      spillover_tmp <- estimate_spillover(ff = ff,
                                          detector = detector,
                                          other_detectors = detectors,
                                          meta = SI_tmp)
      SI[id, colnames(spillover_tmp)] <-  spillover_tmp
    }

    if(estimate_spread){
      if(is.null(comp)) { # Make empty identity matrix with only this detector filled out
        comp_tmp <- diag(length(detectors))
        colnames(comp_tmp) <- rownames(comp_tmp) <- detectors
        comp_tmp[detector, ] <- unlist(spillover_tmp[, grep("Comp", colnames(spillover_tmp))])
      } else {
        comp_tmp <- comp
      }

      spread_tmp <- estimate_spread(ff = ff,
                                    detector = detector,
                                    SI = SI_tmp,
                                    comp = comp_tmp)
      SI[id, colnames(spread_tmp)] <-  spread_tmp
    }

    if (return_cells) {
      set.seed(seed)
      cells <- rbind(
        cells,
        sample_cells(ff,
                     detector,
                     meta = data.frame(ID = id)
        )
      )
    }
  }

  close(pb)

  SI <- indicate_optimal_voltages(SI)

  return(list(
    SI = SI,
    cells = cells
  ))
}

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
    flowCore::keyword(ff, paste0("$P", which(colnames(ff) == detector), "V"))
  )
  meta["Cutoff"] <- cutoff
  meta["Pos_count"] <- sum(pos)
  meta["Neg_count"] <- sum(neg)
  meta["MFI_pos"] <- stats::quantile(flowCore::exprs(ff)[pos, detector], 0.50)
  meta["MFI_neg"] <- stats::quantile(flowCore::exprs(ff)[neg, detector], 0.50)
  meta["q95_neg"] <- stats::quantile(flowCore::exprs(ff)[neg, detector], 0.95)
  meta["q05_neg"] <- stats::quantile(flowCore::exprs(ff)[neg, detector], 0.05)
  meta["rSD"] <- (meta["q95_neg"] - meta["q05_neg"]) / 3.29
  meta["SI"] <- (meta["MFI_pos"] - meta["MFI_neg"]) / (2 * meta["rSD"])

  pData <- Biobase::pData(flowCore::parameters(ff))
  limit <- pData[pData$name == detector, "maxRange"]
  meta["Pctg_OutOfRange"] <- sum(flowCore::exprs(ff)[, detector] >= limit) /
    flowCore::nrow(ff)

  return(meta)
}

#' Estimate spillover for one detector
#' @param ff        FlowFrame. Assumed to be pregated.
#' @param detector  Detector to evaluate. Should be a column name of ff.
#' @param other_detectors Detectors to computer spillover into.
#' @param meta      Dataframe of 1 row with relevant meta information,
#'                  as returned by estimate_SI. Should at least contain
#'                  "Cutoff" column.
#' @export
estimate_spillover <- function(ff,
                               detector,
                               other_detectors,
                               meta){

  pos <- flowCore::exprs(ff)[, detector] >= meta[1, "Cutoff"]
  neg <- flowCore::exprs(ff)[, detector] < meta[1, "Cutoff"]

  for(detector2 in other_detectors){

    # Estimate compensation
    meta[paste0("MFI_pos_", detector2)] <-
      mfi_pos_d2 <- quantile(ff@exprs[pos, detector2], 0.50)
    meta[paste0("MFI_neg_", detector2)] <-
      mfi_neg_d2 <- quantile(ff@exprs[neg, detector2], 0.50)
    meta[paste0("Comp_", detector2)] <- (mfi_pos_d2 - mfi_neg_d2) /
      (meta["MFI_pos"] -  meta["MFI_neg"])
  }

  return(meta)
}

#' Extract a compensation matrix
#'
#' @param SI Result of estimate_brightness with estimate_spillover = TRUE
#' @param singles_of_interest Rownames of the singles you want to include
#'
#' @export
extract_spillover <- function(SI, singles_of_interest){
  detectors <- SI[singles_of_interest, "Detector"]
  comp <- SI[singles_of_interest,
             paste0("Comp_", detectors)]
  colnames(comp) <- gsub("Comp_", "", colnames(comp))
  return(comp)
}

estimate_spread <- function(ff,
                            detector,
                            SI,
                            comp){

  pos <- flowCore::exprs(ff)[, detector] >= SI[["Cutoff"]]
  neg <- flowCore::exprs(ff)[, detector] < SI[["Cutoff"]]

  ff_c <- flowCore::compensate(ff, comp)

  detectors <- colnames(comp)
  for(detector2 in detectors){
    d2_q05_neg <- quantile(ff_c@exprs[neg, detector2], 0.05)
    d2_q50_neg <- quantile(ff_c@exprs[neg, detector2], 0.50)
    d2_q84_neg <- quantile(ff_c@exprs[neg, detector2], 0.84)
    d2_q95_neg <- quantile(ff_c@exprs[neg, detector2], 0.95)
    d2_q50_pos <- quantile(ff_c@exprs[pos, detector2], 0.50)
    d2_q84_pos <- quantile(ff_c@exprs[pos, detector2], 0.84)
    d2_rsd_neg <- (d2_q95_neg - d2_q05_neg) / 3.29
    d2_sigma2_neg <- (d2_q84_neg - d2_q50_neg)^2
    d2_sigma2_pos <- (d2_q84_pos - d2_q50_pos)^2

    if(detector2 != detector & d2_sigma2_pos  > d2_sigma2_neg){
      SI[paste0("Spread_", detector2)] <-
        sqrt(d2_sigma2_pos - d2_sigma2_neg) /
        sqrt(SI["MFI_pos"] -  SI["MFI_neg"])
    } else {
      SI[paste0("Spread_", detector2)] <- 0
    }

    SI[paste0("MFI_pos_c_", detector2)] <- d2_q50_pos
    SI[paste0("MFI_neg_c_", detector2)] <- d2_q50_neg
    if(detector != detector2){
      SI[paste0("SSI_", detector2)] <- (d2_q50_pos - d2_q50_neg) / (2 * d2_rsd_neg)
    } else {
      SI[paste0("SSI_", detector2)] <- 0
    }
  }

  return(SI)
}
