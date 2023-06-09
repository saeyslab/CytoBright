% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_cutoff.R
\name{find_cutoff_flowDensity}
\alias{find_cutoff_flowDensity}
\title{Decide positivity cutoff with flowDensity}
\usage{
find_cutoff_flowDensity(ff, detector, transform = TRUE, ...)
}
\arguments{
\item{ff}{FlowFrame to use}

\item{detector}{Column name of the flow frame}

\item{transform}{Can either be a logical value or a transformList.
If FALSE, no transform is applied. If TRUE (default),
flowcore::estimateLogicle is called, and if this fails, the
default logicleTransform() is applied. If a transformList,
this transformList is applied.}

\item{...}{other parameters to be passed to flowDensity}
}
\description{
Uses the deGate function from the flowDensity package.
Will apply a transformation with estimateLogicle and estimate in the
transformed space, but returns the cutoff in the original space.
}
