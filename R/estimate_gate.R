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
      flowDensity::deGate(ff,
        m,
        all.cuts = TRUE,
        tinypeak.removal = tinypeak.removal
      ),
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
  }


  if (plot) {
    p <- list()
    for (i in seq(1, length(marker_values), 2)) {
      m1 <- names(marker_values)[i]
      if (i + 1 <= length(marker_values)) {
        m2 <- names(marker_values)[i + 1]
      } else {
        m1 <- colnames(ff)[1]
        m2 <- names(marker_values)[i]
      }

      if (!"Original_ID" %in% colnames(ff)) {
        ff <- flowCore::fr_append_cols(
          ff,
          matrix(seq_len(nrow(ff)),
            ncol = 1,
            dimnames = list(NULL, "Original_ID")
          )
        )
      }
      p[[length(p) + 1]] <- filter_plot(ff, ff[selection, ], "", m1, m2) + # ggcyto::autoplot(ff, m1, m2, bins = 128) +
        ggplot2::geom_vline(xintercept = marker_values[[m1]]["min"], color = c("cyan")) +
        ggplot2::geom_vline(xintercept = marker_values[[m1]]["max"], color = c("red")) +
        ggplot2::geom_hline(yintercept = marker_values[[m2]]["min"], color = c("cyan")) +
        ggplot2::geom_hline(yintercept = marker_values[[m2]]["max"], color = c("red")) +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle(plot_title)
    }
    return(list(
      selection = selection,
      plot = p
    ))
  } else {
    return(selection)
  }
}
