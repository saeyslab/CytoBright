% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_plot.R
\name{filter_plot}
\alias{filter_plot}
\title{Plot preprocessing filtering}
\usage{
filter_plot(
  ff_pre,
  ff_post,
  title,
  channel_x,
  channel_y,
  n = 10000,
  id_column = "Original_ID"
)
}
\arguments{
\item{ff_pre}{FlowFrame before filtering}

\item{ff_post}{FlowFrame after filtering}

\item{title}{Plot title}

\item{channel_x}{Channel to display on x-axis}

\item{channel_y}{Channel to display on y-axis}

\item{n}{Maximum number of cells/dots to show}

\item{id_column}{Column to match cells between ff_pre and ff_post}
}
\description{
Plot preprocessing filtering
}
