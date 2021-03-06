% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homogenizationfunctions.R
\name{buishand_range_test}
\alias{buishand_range_test}
\title{Buishand range test}
\usage{
buishand_range_test(
  data_series,
  dp = 2,
  all = FALSE,
  plot = FALSE,
  na.rm = FALSE,
  ...
)
}
\arguments{
\item{data_series}{A climatological time series.}

\item{dp}{Integer indicating the number of decimal places (round) to be used. By default "dp=2".}

\item{all}{A logical value indicating whether to return all Buishand range test values otherwise the absolute maximum is returned.}

\item{plot}{A logical value to return plot of Buishand range test values.}

\item{na.rm}{A logical value indicating whether missing values should be stripped before the computation proceeds.}

\item{...}{Additional arguments to be passed into methods.}
}
\value{
The Buishand range test results.
}
\description{
Buishand range test
}
\author{
Shadrack Kibet
}
