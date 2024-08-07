#' For pregating of a flowFrame
#'
#' @param ff flowFrame to gate
#' @param marker_values List with a vector for every marker to be gated,
#'                      which can contain an estimated minimum threshold (min)
#'                      and the acceptable range to flowDensity cutoff can be in
#'                      (from min - range_min to min + range_min). Same for max
#'                      and max_range. If no possible cut is found in this range,
#'                      the default value is returned (as passed in min or max).
#' @param plot Boolean whether to plot
#' @param plot_title Title for the plot. Default = "Gate"
#' @param tinypeak.removal Passed to flowDensity::deGate, default 1/25
#'
#' @importFrom flowDensity deGate
#' @importFrom flowCore fr_append_cols
#' @export
estimate_gate <- function(ff,
                          marker_values = list(
                            "FSC-A" = c(
                              min = 60000,
                              max = 125000,
                              range_min = 50000,
                              range_max = 50000
                            ),
                            "SSC-A" = c(
                              max = 80000,
                              range_max = 50000
                            )
                          ),
                          plot = FALSE,
                          plot_title = "Gate",
                          tinypeak.removal = 1 / 25) {
  selection <- rep(TRUE, nrow(ff))
  for (m in names(marker_values)) {
    cutoffs <- c(
      tryCatch({flowDensity::deGate(ff,
                                    m,
                                    all.cuts = TRUE,
                                    tinypeak.removal = tinypeak.removal)},
               error = function(e){warning(e); return(NA)}),
      flowDensity::deGate(ff,
                          m,
                          use.upper = TRUE,
                          upper = TRUE
      )
    )
    for (type in c("min", "max")) {
      if (!is.na(marker_values[[m]][type])) {
        default <- marker_values[[m]][type]
        cutoff <- cutoffs[which.min(abs(cutoffs - default))]
        if (!is.na(marker_values[[m]][paste0("range_", type)]) &
          abs(cutoff - default) > marker_values[[m]][paste0("range_", type)]) {
          cutoff <- default
        }
        if (type == "min") selection <- selection & ff@exprs[, m] > cutoff
        if (type == "max") selection <- selection & ff@exprs[, m] < cutoff
        marker_values[[m]][type] <- cutoff
      }
    }
    if(!"plot_min" %in% names(marker_values[[m]])) marker_values[[m]]["plot_min"] <- NA
    if(!"plot_max" %in% names(marker_values[[m]])) marker_values[[m]]["plot_max"] <- NA
  }


  if (plot) {
    p <- list()
    for (i in seq(1, length(marker_values), 2)) {
      m1 <- names(marker_values)[i]
      if (i + 1 <= length(marker_values)) {
        m2 <- names(marker_values)[i + 1]
      } else {
        m1 <- colnames(ff)[1]
        marker_values[[m1]] <- c()
        marker_values[[m1]]["plot_min"] <- NA
        marker_values[[m1]]["plot_max"] <- NA
        m2 <- names(marker_values)[i]
      }


      ff <- flowCore::fr_append_cols(
        ff,
        matrix(seq_len(nrow(ff)),
               ncol = 1,
               dimnames = list(NULL, "Pregating_ID")
        )
      )
      p[[length(p) + 1]] <- filter_plot(ff, ff[selection, ], "", m1, m2, id_column = "Pregating_ID") +
        ggplot2::geom_vline(xintercept = marker_values[[m1]]["min"], color = c("cyan")) +
        ggplot2::geom_vline(xintercept = marker_values[[m1]]["max"], color = c("red")) +
        ggplot2::geom_hline(yintercept = marker_values[[m2]]["min"], color = c("cyan")) +
        ggplot2::geom_hline(yintercept = marker_values[[m2]]["max"], color = c("red")) +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle(plot_title)
      if(!is.na(marker_values[[m1]]["plot_min"]) |
         !is.na(marker_values[[m1]]["plot_max"]))
        p[[length(p)]] <- p[[length(p)]] + ggplot2::xlim(marker_values[[m1]]["plot_min"],
                                                         marker_values[[m1]]["plot_max"])
      if(!is.na(marker_values[[m2]]["plot_min"]) |
         !is.na(marker_values[[m2]]["plot_max"]))
        p[[length(p)]] <- p[[length(p)]] + ggplot2::ylim(marker_values[[m2]]["plot_min"],
                                                         marker_values[[m2]]["plot_max"])
    }
    return(list(
      selection = selection,
      plot = p
    ))
  } else {
    return(selection)
  }
}
