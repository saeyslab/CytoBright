#' Find plateau from numeric vector
#'
#' Takes every point as possible split,
#' fits linear model on left side and horizontal line on right side,
#' and looks for the point with minimal residuals.
#'
#' @param data Vector with values
#'
#' @returns index of the start of the plateau
#'
#' @importFrom stats lm median
#'
#' @export
find_plateau <- function(data) {
  data <- as.data.frame(cbind(seq_along(data), data))
  colnames(data) <- c("X", "Y")

  data <- data[!is.na(data$Y), ]

  n <- nrow(data)
  if (n > 1) {
    min_r <- Inf
    optimal <- 1
    for (i in 2:(n - 1)) {
      f1 <- stats::lm(Y ~ X, data[1:i, ])
      f2 <- stats::median(data[i:n, "Y"])
      f2_residuals <- f2 - data[i:n, "Y"]
      r <- sum(abs(c(f1$residuals, f2_residuals)))
      if (r < min_r) {
        min_r <- r
        optimal <- i
      }
    }

    return(data$X[optimal])
  } else {
    return(1)
  }
}
