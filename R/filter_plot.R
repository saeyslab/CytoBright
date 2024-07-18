#' Plot preprocessing filtering
#'
#' @param ff_pre  FlowFrame before filtering
#' @param ff_post FlowFrame after filtering
#' @param title   Plot title
#' @param channel_x Channel to display on x-axis
#' @param channel_y Channel to display on y-axis
#' @param n         Maximum number of cells/dots to show
#' @param id_column Column to match cells between ff_pre and ff_post
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
                        n = 10000,
                        id_column = "Original_ID") {

  if (!id_column %in% colnames(flowCore::exprs(ff_pre))) {
    ff_pre <- flowCore::fr_append_cols(ff_pre,
                                       matrix(seq_len(nrow(ff_pre)),
                                              ncol = 1,
                                              dimnames = list(NULL,
                                                              c(id_column))))
  }
  df <- data.frame(
    x = flowCore::exprs(ff_pre)[, channel_x],
    y = flowCore::exprs(ff_pre)[, channel_y],
    selected = flowCore::exprs(ff_pre)[, id_column] %in%
      flowCore::exprs(ff_post)[, id_column]
  )

  i <- sample(nrow(df), min(n, nrow(df)))
  p <- ggplot2::ggplot(df[i, ],
                       ggplot2::aes(x = .data$x,
                                    y = .data$y,
                                    color = .data$selected)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
    ggplot2::xlab(paste(flowCore::getChannelMarker(ff_pre, channel_x), collapse = " ")) +
    ggplot2::ylab(paste(flowCore::getChannelMarker(ff_pre, channel_y), collapse = " ")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle(paste0(round(flowCore::nrow(ff_post) / flowCore::nrow(ff_pre) * 100, 2), "% ", title))
  return(p)
}
