#' Estimate SI from single stains
#'
#' @param single_stains Dataframe with at least the following columns:
#'                      ID, Fluorochrome, Detector, File.
#'                      ID is assumed to be unique.
#'                      Detector should correspond with the colnames of the
#'                      flow frames.
#'                      File should be full path to the single stain fcs files.
#' @param preprocessing_parameters Parameters for pregating done by flowDensity.
#'                      Named list with booleans (compensate, removeMargins,
#'                      removeDoublets) and "pregate", which is a list itself,
#'                      and contains per gate a list marker_values.
#'                      The name corresponding to the channel.
#'
#' @param return_cells  Boolean. If true, a list with cell values are returned
#' @param seed          Seed for reproducability
#' @param comp          Compensation matrix.
#' @param transform Can either be a logical value or a transformList.
#'                  If FALSE, no transform is applied. If TRUE (default),
#'                  flowcore::estimateLogicle is called, and if this fails, the
#'                  default logicleTransform() is applied. If a transformList,
#'                  this transformList is applied. Note that the value is
#'                  returned in the original space.
#' @param unstained     If this parameter is provided, it should contain a path
#'                      to an fcs file. The unstained population will then be
#'                      taken from this file rather than the single stain.
#'
#' @importFrom flowCore read.FCS transform transformList arcsinhTransform colnames compensate
#' @importFrom flowDensity deGate
#' @importFrom PeacoQC RemoveMargins RemoveDoublets
#' @importFrom FlowSOM AggregateFlowFrames
#'
#' @export
estimate_brightness <- function(single_stains,
                                preprocessing_parameters =
                                  list(
                                    compensate = FALSE,
                                    removeMargins = TRUE,
                                    removeDoublets = TRUE,
                                    pregate = list(list(
                                      marker_values = list(
                                        "FSC-A" = c(
                                          min = 60000,
                                          range_min = 50000,
                                          max = 125000,
                                          range_max = 50000
                                        ),
                                        "SSC-A" = c(
                                          max = 80000,
                                          range_max = 50000
                                        )
                                      )
                                    )),
                                    pregate_tf = FALSE
                                  ),
                                return_cells = TRUE,
                                seed = 1,
                                comp = NULL,
                                transform = TRUE,
                                unstained = NULL,
                                estimate_spillover = FALSE) {
  rownames(single_stains) <- single_stains$ID
  if (is.null(single_stains$Group)) single_stains$Group <- single_stains$ID
  groups <- unique(single_stains$Group)
  detectors <- unique(single_stains$Detector)

  pb <- utils::txtProgressBar(min = 0, max = length(groups), style = 3)

  values_of_interest <- c(
    "Cutoff", "MFI_pos", "MFI_neg",
    "Min_neg", "Max_neg", "rSD",
    "SI", "Voltage"
  )

  if(estimate_spillover){
    values_for_comp <- c(paste0("MFI_pos_", detectors),
                         paste0("MFI_neg_", detectors),
                         paste0("Comp_", detectors),
                         paste0("Spread_", detectors))
  } else {
    values_for_comp <- c()
  }

  SI <- data.frame(
    matrix(
      nrow = 0,
      ncol = ncol(single_stains) +
        length(values_of_interest) +
        length(values_for_comp),
      dimnames = list(
        NULL,
        c(
          colnames(single_stains),
          values_of_interest,
          values_for_comp
        )
      )
    ),
    check.names = FALSE
  )

  pregating_plots <- list()

  cells <- data.frame(matrix(NA,
    nrow = 0,
    ncol = 2,
    dimnames = list(NULL, c("ID", "Value"))
  ))

  if (!is.null(unstained)) {
    ff_unstained <- do.call(
      prep_FCS,
      c(
        list(file = unstained),
        preprocessing_parameters
      )
    )
  }

  for (group in groups) {
    utils::setTxtProgressBar(pb, which(groups == group))

    row_ids <- which(single_stains$Group == group)

    fluor <- single_stains[row_ids[1], "Fluorochrome"]
    detector <- single_stains[row_ids[1], "Detector"]

    if (length(row_ids) > 1) {
      files <- single_stains[row_ids, "File"]
      set.seed(seed)
      ff <- FlowSOM::AggregateFlowFrames(files,
        cTotal = 3000000,
        truncate_max_range = FALSE,
        silent = TRUE
      )
    } else {
      file <- single_stains[row_ids, "File"]
      ff <- flowCore::read.FCS(file, truncate_max_range = FALSE)
    }


    ff_prep <- do.call(
      prep_FCS,
      c(
        list(file = ff),
        preprocessing_parameters
      )
    )
    ff <- ff_prep$flowFrame
    pregating_plots[[group]] <- ff_prep$plot

    if ("File" %in% colnames(ff)) {
      subsets <- ff@exprs[, "File"]
    } else {
      subsets <- rep(1, nrow(ff))
    }

    for (subset in unique(subsets)) {
      selection <- subsets == subset
      if (length(unique(subsets)) > 1) {
        if (length(row_ids) == length(unique(subsets))) {
          sub_id <- single_stains[row_ids[subset], "ID"]
          SI[sub_id, colnames(single_stains)] <- single_stains[sub_id, ]
        } else {
          sub_id <- paste0(group, "_", subset)
          SI[sub_id, colnames(single_stains)] <- single_stains[group, ]
          SI[sub_id, "ID"] <- sub_id
        }
      } else {
        sub_id <- group
        SI[sub_id, colnames(single_stains)] <- single_stains[sub_id, ]
      }

      ff_tmp <- ff[selection, ]
      if(length(row_ids) > 1){
        # restore voltage which was lost in concatenation
        # Thanks to Juan Hernandez
        voltage_keyword <- paste0("$P", which(colnames(ff) == detector), "V")
        orig_file <- flowCore::read.FCS(single_stains[row_ids[subset], "File"],
                                        which.lines = 1,
                                        truncate_max_range = FALSE)
        flowCore::keyword(ff_tmp)[[voltage_keyword]] <-
          flowCore::keyword(orig_file)[[voltage_keyword]]
      }

      if (!is.null(unstained)) {
        cutoff <- find_cutoff_FMO(ff_unstained$flowFrame, detector, 0.995)

        flowCore::exprs(ff_tmp) <- rbind(
          flowCore::exprs(ff_tmp),
          flowCore::exprs(ff_unstained$flowFrame)
        )
      } else {
        cutoff <- find_cutoff_flowDensity(ff_tmp,
                                          detector = detector,
                                          transform = transform)
      }

      SI_tmp <- estimate_SI(
        ff = ff_tmp,
        detector = detector,
        cutoff = cutoff
      )

      SI[sub_id, colnames(SI_tmp)] <- SI_tmp


      if(estimate_spillover){
        spillover_tmp <- estimate_spillover(ff = ff_tmp,
                                            detector = detector,
                                            other_detectors = detectors,
                                            SI = SI_tmp)
        SI[sub_id, colnames(spillover_tmp)] <-  spillover_tmp

        if(is.null(comp)) { # Make empty identity matrix with only this detector filled out
          comp_tmp <- diag(length(detectors))
          colnames(comp_tmp) <- rownames(comp_tmp) <- detectors
          comp_tmp[detector, ] <- unlist(spillover_tmp[,grep("Comp", colnames(spillover_tmp))])
        } else {
          comp_tmp <- comp

          spread_tmp <- estimate_spread(ff = ff_tmp,
                                        detector = detector,
                                        SI = SI_tmp,
                                        comp = comp_tmp)
          SI[sub_id, colnames(spread_tmp)] <-  spread_tmp
        }
      }

      if (return_cells) {
        set.seed(seed)
        cells <- rbind(
          cells,
          sample_cells(ff_tmp,
            detector,
            meta = data.frame(ID = sub_id)
          )
        )
      }
    }
  }
  close(pb)

  SI <- indicate_optimal_voltages(SI)

  return(list(
    SI = SI,
    cells = cells,
    pregating_plots = pregating_plots
  ))
}

#' Estimate spillover for one detector
#' @param ff        FlowFrame. Assumed to be pregated.
#' @param detector  Detector to evaluate. Should be a column name of ff.
#' @param other_detectors Detectors to computer spillover into.
#' @param SI        Dataframe of 1 row with relevant meta information,
#'                  as returned by estimate_SI. Should at least contain
#'                  "Cutoff" column.
#' @export
estimate_spillover <- function(ff,
                               detector,
                               other_detectors,
                               SI){

  pos <- flowCore::exprs(ff)[, detector] >= SI[["Cutoff"]]
  neg <- flowCore::exprs(ff)[, detector] < SI[["Cutoff"]]

  for(detector2 in other_detectors){
    # Estimate compensation
    SI[paste0("MFI_pos_", detector2)] <-
      mfi_pos_d2 <- quantile(ff@exprs[pos, detector2], 0.50)
    SI[paste0("MFI_neg_", detector2)] <-
      mfi_neg_d2 <- quantile(ff@exprs[neg, detector2], 0.50)
    SI[paste0("Comp_", detector2)] <- (mfi_pos_d2 - mfi_neg_d2) /
      (SI["MFI_pos"] -  SI["MFI_neg"])
  }

  return(SI)
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
    d2_q50_neg <- quantile(ff_c@exprs[neg, detector2], 0.50)
    d2_q84_neg <- quantile(ff_c@exprs[neg, detector2], 0.84)
    d2_q50_pos <- quantile(ff_c@exprs[pos, detector2], 0.50)
    d2_q84_pos <- quantile(ff_c@exprs[pos, detector2], 0.84)
    d2_sigma2_neg      <- (d2_q84_neg - d2_q50_neg)^2
    d2_sigma2_pos  <- (d2_q84_pos - d2_q50_pos)^2

    if(detector2 != detector & d2_sigma2_pos  > d2_sigma2_neg){
      SI[paste0("Spread_", detector2)] <-
        sqrt(d2_sigma2_pos - d2_sigma2_neg) /
        sqrt(SI["MFI_pos"] -  SI["MFI_neg"])
    } else {
      SI[paste0("Spread_", detector2)] <- 0
    }
  }

  return(SI)
}
