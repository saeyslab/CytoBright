#' Sample cells in long format
#'
#' @param ff flowFrame
#' @param detector Column name of flowFrame
#' @param meta Dataframe with extra metadata describing the flowFrame. Will be
#'             added as additional columns to the result
#' @param n Number of cells to sample at most. If flowFrame contains fewer cells,
#'          nrow(ff) will be sampled.
#'
#'
#' @importFrom flowCore exprs
sample_cells <- function(ff, detector, meta, n = 6000) {
  cells_tmp <- data.frame(
    Detector = unname(detector),
    Value = sample(
      flowCore::exprs(ff)[, detector],
      min(n, nrow(ff))
    )
  )
  cells <- cbind(meta, cells_tmp)
  return(cells)
}
