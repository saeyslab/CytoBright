#' Plot preprocessing filtering
#'
#' @param ff_pre  FlowFrame before filtering
#' @param ff_post FlowFrame after filtering
#' @param title   Plot title
#' @param channel_x Channel to display on x-axis
#' @param channel_y Channel to display on y-axis
#' @param n         Maximum number of cells/dots to show
#'
#' @importFrom flowCore exprs fr_append_cols getChannelMarker nrow
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes xlab ylab theme_minimal theme ggtitle
#' @export
filter_plot <- function(ff_pre,
                        ff_post,
                        title,
                        channel_x,
                        channel_y,
                        n = 10000) {
  df <- data.frame(
    x = flowCore::exprs(ff_pre)[, channel_x],
    y = flowCore::exprs(ff_pre)[, channel_y]
  )
  i <- sample(nrow(df), min(n, nrow(df)))
  if (!"Original_ID" %in% colnames(flowCore::exprs(ff_pre))) {
    ff_pre <- flowCore::fr_append_cols(
      ff_pre,
      matrix(seq_len(nrow(df)),
        ncol = 1,
        dimnames = list(NULL, c("Original_ID"))
      )
    )
  }
  p <- ggplot2::ggplot(df[i, ], ggplot2::aes(
    x = .data$x,
    y = .data$y
  )) +
    ggplot2::geom_point(
      size = 0.5,
      color = ifelse(flowCore::exprs(ff_pre)[i, "Original_ID"] %in%
        flowCore::exprs(ff_post)[, "Original_ID"], "blue", "red")
    ) +
    ggplot2::xlab(paste(flowCore::getChannelMarker(ff_pre, channel_x), collapse = " ")) +
    ggplot2::ylab(paste(flowCore::getChannelMarker(ff_pre, channel_y), collapse = " ")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle(paste0(round(flowCore::nrow(ff_post) / flowCore::nrow(ff_pre) * 100, 2), "% ", title))
  return(p)
}
