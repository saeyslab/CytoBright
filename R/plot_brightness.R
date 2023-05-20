#' Plot SI as computed by estimate_brightness
#'
#' @param SI_optimal SI result to plot
#' @param cells      As returned by estimate_brightness, to show signal intensity
#' @param cofactor   asinh cofactor used to plot signal intensity
#' @param order      Column of SI_optimal used for sorting
#' @param decreasing Argument for sorting
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_violin geom_text
#'                     element_text
#'                     theme theme_minimal theme_void
#'                     ggtitle xlab ylab
#'                     scale_x_discrete scale_fill_manual
#' @importFrom patchwork plot_layout
#' @export
plot_brightness <- function(SI_optimal,
                            cells = NULL,
                            cofactor = 150,
                            order = "SI",
                            decreasing = TRUE) {
  if (!is.null(order)) {
    SI_optimal$ID <- factor(SI_optimal$ID,
      levels = SI_optimal$ID[order(SI_optimal[, order],
        decreasing = !decreasing
      )]
    )
    # reverse decreasing because we look top to bottom but ggplot bottom-to-top
  }

  if (is.null(SI_optimal$BrightnessLevel)) SI_optimal$BrightnessLevel <- NA
  p_SI <- ggplot2::ggplot(SI_optimal) +
    ggplot2::geom_col(ggplot2::aes(
      y = .data$ID,
      x = .data$SI,
      fill = as.character(.data$BrightnessLevel)
    )) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = .data$ID,
        x = .data$SI + 3, label = round(.data$SI)
      ),
      size = 2
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text( # angle = 90, hjust = 0.5,vjust=1,
      face = "bold"
    )) +
    ggplot2::scale_fill_manual(
      values = c(
        "1" = "#2166ac",
        "2" = "#67a9cf",
        "3" = "#fddbc7",
        "4" = "#ef8a62",
        "5" = "#b2182b"
      ),
      na.value = "lightgrey",
      name = "Theoretical brightness"
    ) +
    ggplot2::ylab("") + # Fluorochrome
    ggplot2::xlab("Stain index") +
    ggplot2::scale_y_discrete(
      breaks = SI_optimal$ID,
      labels = SI_optimal$ID,
      position = "right"
    ) +
    ggplot2::scale_x_continuous(position = "top")
  plots <- list("SI" = p_SI)
  heights <- c(3)

  if (!is.null(cells)) {
    cells_subset <- cells %>% dplyr::filter(.data$ID %in% rownames(SI_optimal))
    cells_subset$ID <- factor(cells_subset$ID, levels = levels(SI_optimal$ID))
    p_signal <- ggplot2::ggplot(data.frame(
      ID = cells_subset[, "ID"],
      Value = cells_subset[, "Value"]
    )) +
      ggplot2::geom_violin(
        ggplot2::aes(
          x = asinh(.data$Value / cofactor),
          y = .data$ID,
          group = .data$ID
        ),
        fill = "black",
        size = 0.5
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          y = .data$ID,
          x = asinh(.data$MFI_pos / cofactor)
        ),
        data = SI_optimal,
        col = "red"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          y = .data$ID,
          x = asinh(.data$MFI_neg / cofactor)
        ),
        data = SI_optimal,
        col = "cyan"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          y = .data$ID,
          x = asinh(.data$Max_neg / cofactor)
        ),
        data = SI_optimal,
        col = "cyan", size = 1
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          y = .data$ID,
          x = asinh(.data$Min_neg / cofactor)
        ),
        data = SI_optimal,
        col = "cyan", size = 1
      ) +
      ggplot2::xlim(c(-1, 10)) +
      ggplot2::ylab("") +
      ggplot2::xlab(paste0("asinh(x/", cofactor, ")")) +
      ggplot2::scale_x_continuous(position = "top")


    heights <- c(heights, 3)
    plots[["signal"]] <- p_signal
  }

  if (!is.null(SI_optimal$Voltage)) {
    if (!is.null(SI_optimal$Ref_Voltage)) {
      p_voltages <- ggplot2::ggplot(SI_optimal) +
        ggplot2::geom_point(ggplot2::aes(y = .data$ID, x = as.numeric(.data$Voltage))) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::xlab("Voltage") +
        ggplot2::ylab("")

      p_voltages <- p_voltages +
        ggplot2::geom_point(
          ggplot2::aes(
            x = as.numeric(.data$RefVoltage),
            y = .data$ID
          ),
          col = "grey"
        )
    } else {
      p_voltages <- ggplot2::ggplot(SI_optimal) +
        ggplot2::geom_text(ggplot2::aes(y = .data$ID, x = 1, label = .data$Voltage),
          size = 3,
        ) + # angle = 90, hjust = 0.5, vjust = 0
        ggplot2::theme_void() +
        ggplot2::xlab("Voltage") +
        ggplot2::theme(axis.title.x = element_text()) +
        ggplot2::scale_x_continuous(position = "top")
    }

    heights <- c(heights, 1)
    plots[["voltages"]] <- p_voltages
  }

  if (!is.null(SI_optimal$`Detector`)) {
    plots[["detector"]] <- ggplot2::ggplot(SI_optimal) +
      ggplot2::geom_text(
        ggplot2::aes(
          y = .data$ID, x = 1,
          label = gsub("-A", "", .data$Detector)
        ),
        size = 3
      ) + # angle = 90, vjust = 0.5, hjust = 0,
      ggplot2::theme_void() +
      ggplot2::xlab("Detector") +
      ggplot2::theme(axis.title.x = element_text()) +
      ggplot2::scale_x_continuous(position = "top")

    heights <- c(heights, 2)
  }

  if (!is.null(SI_optimal$`Laser`)) {
    plots[["laser"]] <- ggplot2::ggplot(SI_optimal) +
      ggplot2::geom_text(
        ggplot2::aes(
          y = .data$ID, x = 0,
          label = .data$`Laser`
        ),
        size = 3, fontface = "bold"
      ) +
      ggplot2::geom_rect(aes(
        ymin = as.numeric(.data$ID) - 0.3,
        ymax = as.numeric(.data$ID) + 0.3,
        xmin = -0.05, xmax = 0.1,
        fill = .data$`Laser`
      ), alpha = 0.8) +
      ggplot2::geom_text(
        ggplot2::aes(
          y = .data$ID, x = 0,
          label = .data$`Laser`
        ),
        # angle = 90, hjust = 0.5, vjust = 0,
        size = 3, fontface = "bold", col = "black"
      ) +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")) +
      ggplot2::scale_fill_manual(values = c(
        "R" = "#cc0000",
        "B" = "#0085cc",
        "V" = "#de43c7",
        "YG" = "#32a60f",
        "UV" = "#97889e"
      )) +
      # ggplot2::theme(legend.position = "none") +
      ggplot2::guides(fill = "none") +
      ggplot2::theme_void() +
      ggplot2::xlab("Laser") +
      ggplot2::theme(axis.title.x = element_text()) +
      ggplot2::scale_x_continuous(position = "top")

    heights <- c(heights, 1)
  }


  if (!is.null(SI_optimal$Clone)) {
    plots[["clone"]] <- ggplot2::ggplot(SI_optimal) +
      ggplot2::geom_text(ggplot2::aes(y = .data$ID, x = 1, label = .data$Clone),
        size = 3
      ) + # angle = 90, hjust = 0.5, vjust = 0,
      ggplot2::theme_void() +
      ggplot2::xlab("Clone") +
      ggplot2::theme(axis.title.x = element_text()) +
      ggplot2::scale_x_continuous(position = "top")

    heights <- c(heights, 1)
  }

  if (!is.null(SI_optimal$Brand)) {
    plots[["brands"]] <- ggplot2::ggplot(SI_optimal) +
      ggplot2::geom_text(ggplot2::aes(y = .data$ID, x = 1, label = .data$Brand),
        size = 3
      ) + # angle = 90, hjust = 0.5, vjust = 0,
      ggplot2::theme_void() +
      ggplot2::xlab("Brand") +
      ggplot2::theme(axis.title.x = element_text()) +
      ggplot2::scale_x_continuous(position = "top")

    heights <- c(heights, 1)
  }

  plots[["top_axis"]] <- ggplot2::ggplot(SI_optimal) +
    ggplot2::geom_text(ggplot2::aes(y = .data$ID, x = 1, label = .data$ID),
      size = 3, fontface = "bold"
    ) + # hjust = 0,angle = 90, hjust = 0.5,
    # ggplot2::theme(plot.margin = ggplot2::margin(0,0,0,0, "pt")) +
    ggplot2::theme_void() +
    ggplot2::xlab("Fluor") +
    ggplot2::theme(axis.title.x = element_text()) +
    ggplot2::scale_x_continuous(position = "top")
  heights <- c(heights, 2)

  patchwork::wrap_plots(rev(plots), nrow = 1) +
    patchwork::plot_layout(widths = rev(heights), guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}
