#' Preprocess FCS file
#'
#' @param file Path to fcs file or flowFrame.
#'             Will be read with truncate_max_range = FALSE. Additional parameters
#'             for read.FCS will be passed along.
#' @param compensate  Can be FALSE (no compensation), a character (keyword in
#'                    the flowFrame containing the compensation matrix) or a
#'                    numeric matrix (compensation matrix to use).
#'                    Default = "SPILL"
#' @param removeMargins  Can be FALSE (no margin removal), TRUE (margin removal
#'                       for all channels except those containing File or
#'                       Original_ID) or a character vector (containing the
#'                       names of the channels to remove the margins for).
#'                       Default = TRUE
#' @param removeDoublets Boolean, whether or not to remove
#' @param pregate Pregating information
#' @param pregate_tf TransformList to use in pregating.
#' @param ... Extra arguments to pass to read.FCS
#'
#' @importFrom methods is
#' @importFrom flowCore read.FCS keyword compensate getChannelMarker
#'                      estimateLogicle transformList transform
#' @importFrom PeacoQC RemoveMargins RemoveDoublets
#' @importFrom utils capture.output
#' @export
prep_FCS <- function(file,
                     compensate = "SPILL",
                     removeMargins = TRUE,
                     removeDoublets = TRUE,
                     pregate = TRUE,
                     pregate_tf = TRUE,
                     ...) {
  if (is.character(file)) {
    ff <- flowCore::read.FCS(file,
      truncate_max_range = FALSE,
      ...
    )
  } else if (methods::is(file, "flowFrame")) {
    ff <- file
  } else {
    stop("File should be a path to an fcs file or a flowFrame")
  }

  to_plot <- list()

  if (!isFALSE(removeMargins)) {
    if (isTRUE(removeMargins)) {
      removeMargins <- grep("File|Original_ID",
        colnames(ff),
        invert = TRUE, value = TRUE
      )
    }
    ff_m <- PeacoQC::RemoveMargins(ff, removeMargins)
    to_plot[[length(to_plot) + 1]] <- list(
      ff_pre = ff,
      ff_post = ff_m,
      title = "Non-margin events",
      channels = c("FSC-A", "SSC-A")
    )
  } else {
    ff_m <- ff
  }

  if (isTRUE(removeDoublets)) {
    ff_s <- PeacoQC::RemoveDoublets(ff_m)
    to_plot[[length(to_plot) + 1]] <- list(
      ff_pre = ff_m,
      ff_post = ff_s,
      title = "Singlets",
      channels = c("FSC-A", "FSC-H")
    )
  } else {
    ff_s <- ff_m
  }


  plot_list <- list()
  for (plot in to_plot) {
    plot_list[[length(plot_list) + 1]] <- filter_plot(
      ff_pre = plot$ff_pre,
      ff_post = plot$ff_post,
      title = plot$title,
      channel_x = plot$channels[1],
      channel_y = plot$channels[2]
    )
  }

  if (!isFALSE(compensate)) {
    if (is.character(compensate)) {
      compensate <- flowCore::keyword(ff, compensate)[[1]]
    }
    ff_s <- flowCore::compensate(ff_s, compensate)
  }

  if (!isFALSE(pregate)) {
    if (isTRUE(pregate)) {
      marker_values_1 <- list(marker_values = list(
        "FSC-A" = c(
          min = 30000,
          range_min = 10000,
          max = 100000,
          range_max = 10000
        ),
        "SSC-A" = c(
          max = 70000,
          range_max = 10000
        )
      ))

      LD_channel <- flowCore::getChannelMarker(ff, "LD")[, "name"]
      marker_values_2 <- list(marker_values = list())
      marker_values_2[["marker_values"]][[LD_channel]] <- c(
        max = 2,
        range_max = 0.5
      )

      pregate <- list(marker_values_1, marker_values_2)
    }

    if (!isFALSE(pregate_tf)) {
      if (isTRUE(pregate_tf)) {
        pregate_tf <- tryCatch(
          {
            to_transform <- grep("FSC|SSC",
              sapply(pregate, function(x) names(x$marker_values)),
              invert = TRUE, value = TRUE
            )
            flowCore::estimateLogicle(ff_s, to_transform)
          },
          error = function(e) {
            warning(
              "Default logicle parameters for ",
              paste(to_transform, collapse = ",")
            )
            flowCore::transformList(
              to_transform,
              flowCore::logicleTransform()
            )
          }
        )
      }
      ff_t <- flowCore::transform(ff_s, pregate_tf)
    } else {
      ff_t <- ff_s
    }

    ff_gated <- ff_s
    gate_res <- list()
    for (gate in pregate) {
      n <- length(gate_res) + 1
      o <- utils::capture.output(
        gate_res[[n]] <- do.call(
          estimate_gate,
          c(
            list(ff_t),
            gate,
            list(plot = TRUE)
          )
        )
      )
      selection <- gate_res[[n]]$selection
      gate_res[[n]]$plot <- gate_res[[n]]$plot[[1]] +
        ggplot2::ggtitle(paste0(
          round(100 * sum(selection) / length(selection), 2),
          "% ",
          paste(names(gate$marker_values), collapse = " - ")
        ))
      ff_t <- ff_t[selection, ]
      ff_gated <- ff_gated[selection, ]
    }

    plot_list <- c(plot_list, lapply(gate_res, function(x) x$plot))
  } else {
    ff_gated <- ff_s
  }

  return(list(
    flowFrame = ff_gated,
    plot = plot_list,
    tf = pregate_tf
  ))
}
