#' Decide positivity cutoff based on quantile
#'
#' Mainly useful when using an FMO or unstained sample.
#'
#' @param ff FlowFrame to use
#' @param detector Column name of the flow frame
#' @param quantile Quantile used for the cutoff. Default = 0.9999
#'
#' @importFrom flowCore exprs
#'
#' @export
find_cutoff_FMO <- function(ff,
                            detector,
                            quantile = 0.9999) {
  cutoff <- quantile(flowCore::exprs(ff)[, detector], quantile)
  return(unname(cutoff))
}

#' Decide positivity cutoff with flowDensity
#'
#' Uses the deGate function from the flowDensity package.
#' Will apply a transformation with estimateLogicle and estimate in the
#' transformed space, but returns the cutoff in the original space.
#'
#' @param ff FlowFrame to use
#' @param detector Column name of the flow frame
#' @param transform Can either be a logical value or a transformList.
#'                  If FALSE, no transform is applied. If TRUE (default),
#'                  flowcore::estimateLogicle is called, and if this fails, the
#'                  default logicleTransform() is applied. If a transformList,
#'                  this transformList is applied.
#' @param ... other parameters to be passed to flowDensity
#'
#' @importFrom methods is
#' @importFrom flowDensity deGate
#' @importFrom flowCore estimateLogicle transform inverseLogicleTransform
#'
#' @export
find_cutoff_flowDensity <- function(ff, detector, transform = TRUE, ...) {
  if (!transform) {
    ff_t <- ff
  } else {
    if (!methods::is(transform, "transformList")) {
      transform <- tryCatch(
        flowCore::estimateLogicle(ff, detector, type = "data"),
        error = function(e) {
          warning(
            "Default logicle parameters for ",
            paste(detector, collapse = ",")
          )
          flowCore::transformList(
            detector,
            flowCore::logicleTransform()
          )
        }
      )
    }

    ff_t <- flowCore::transform(ff, transform)
  }

  cutoff <- flowDensity::deGate(ff_t, detector,
    tinypeak.removal = 0.2,
    upper = TRUE,
    ...
  )

  if (is(transform, "transformList")) {
    q <- stats::ecdf(flowCore::exprs(ff_t)[, detector])(cutoff)
    cutoff <- stats::quantile(flowCore::exprs(ff)[, detector], q)
  }
  return(cutoff)
}
