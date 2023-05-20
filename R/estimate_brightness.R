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
                                unstained = NULL) {
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

  values_for_comp <- c()
  # values_for_comp <- c(paste0("MFI_pos_", detectors),
  #                      paste0("MFI_neg_", detectors),
  #                      paste0("Comp_", detectors),
  #                      paste0("Spread_", detectors))

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

      if (!is.null(unstained)) {
        cutoff <- find_cutoff_FMO(ff_unstained$flowFrame, detector, 0.995)

        ff_tmp <- ff[selection, ]
        flowCore::exprs(ff_tmp) <- rbind(
          flowCore::exprs(ff_tmp),
          flowCore::exprs(ff_unstained$flowFrame)
        )
      } else {
        cutoff <- find_cutoff_flowDensity(ff[selection, ],
          detector = detector
        )

        ff_tmp <- ff[selection, ]
      }


      SI_tmp <- estimate_SI(
        ff = ff_tmp,
        detector = detector,
        cutoff = cutoff
      )

      SI[sub_id, colnames(SI_tmp)] <- SI_tmp


      # spillover_tmp <- estimate_spillover(ff = ff_tmp,
      #                                     detector = detector,
      #                                     other_detectors = detectors,
      #                                     SI = SI_tmp)
      # SI[sub_id, colnames(spillover_tmp)] <-  spillover_tmp
      #
      # if(is.null(comp)) { # Make empty identity matrix with only this detector filled out
      #   # comp_tmp <- diag(length(detectors))
      #   # colnames(comp_tmp) <- rownames(comp_tmp) <- detectors
      #   # comp_tmp[detector, ] <- unlist(spillover_tmp[,grep("Comp", colnames(spillover_tmp))])
      # } else {
      #   comp_tmp <- comp
      #
      #   spread_tmp <- estimate_spread(ff = ff_tmp,
      #                                 detector = detector,
      #                                 SI = SI_tmp,
      #                                 comp = comp_tmp)
      #   SI[sub_id, colnames(spread_tmp)] <-  spread_tmp
      # }


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

  best_per_group <- tapply(
    SI$SI[SI$Pctg_OutOfRange < 0.01],
    SI$Group[SI$Pctg_OutOfRange < 0.01],
    find_plateau
  )
  best_per_group <- sapply(unique(SI$Group), function(x) {
    SI$ID[which(SI$Group == x)][best_per_group[x]]
  })

  return(list(
    SI = SI,
    cells = cells,
    pregating_plots = pregating_plots,
    best_per_group = best_per_group
  ))
}
