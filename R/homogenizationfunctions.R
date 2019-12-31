
# Homogenization tests ----------------------------------------------------

#' Buishand range test
#'
#' @param na.rm A logical value indicating whether missing values should be stripped before the computation proceeds.
#' @param data_series A climatological time series.
#' @param ... Additional arguments to be passed into methods.
#'
#' @return The Buishand range test results.
#' @author Shadrack Kibet
#' @export
#'
#' @examples
buishand_range_test <- function(data_series, na.rm = FALSE, ...) {
  if(missing(data_series)) stop("argument 'data series' is missing with no default")
  n <- length(data_series)
  s2 <- 1:n
  s <- sapply(s2, function(k) {
    m <- summary_mean(data_series, na.rm = na.rm)
    b <- summary_sum(data_series[1:k] - m, na.rm = na.rm)
    return(b)
  })
  s / summary_sd(data_series, na.rm = na.rm) / sqrt(n)
}

#' Craddock test
#'
#' @param na.rm A logical value indicating whether missing values should be stripped before the computation proceeds.
#' @param dp Integer indicating the number of decimal places (round) to be used. By default "dp=2".
#' @param data_col1 A candidate climatological series.
#' @param data_col2 A homogeneous reference series.
#' @param ... Additional arguments to be passed into methods.
#'
#' @return Craddock test results.
#' @author Shadrack Kibet
#' @export
#'
#' @examples
craddock_test <- function(data_col1, data_col2, na.rm = FALSE, dp = 2, ...) {
  if (missing(data_col1)||missing(data_col2)) stop("both columns 'data_col1' and 'data_col2' must be supplied")
  mean_col1 <- summary_mean(data_col1, na.rm = na.rm)
  mean_col2 <- summary_mean(data_col2, na.rm = na.rm)
  c2C <- data_col2 + mean_col1 - mean_col2
  col_diff <- data_col1 - c2C
  round(cumsum(col_diff), digits = dp)
}

#' Standard Normal Homogeneity test(SNHT)
#'
#' @param na.rm A logical value indicating whether missing values should be stripped before the computation proceeds.
#' @param dp Integer indicating the number of decimal places (round) to be used. By default "dp=2".
#' @param data_series A climatological time series.
#' @param ... Additional arguments to be passed into methods.
#'
#' @return The SNHT test statistic.
#' @author Shadrack Kibet
#' @export
#'
#' @examples
snht <- function(data_series, na.rm = FALSE, dp = 2, ...) {
  if(missing(data_series)) stop("argument 'data series' is missing with no default")
  tsn <- (data_series - summary_mean(data_series, na.rm = na.rm)) / summary_sd(data_series, na.rm = na.rm)
  n <- length(data_series)
  zvec <- vector("numeric", n - 1)
  for (v in 1:(n - 1)) {
    zvec[v] <- v * (summary_mean(tsn[1:v], na.rm = na.rm))^2 +
      (n - v) * (summary_mean(tsn[(v + 1):n], na.rm = na.rm))^2
  }
  round(zvec, digits = dp)
}

#' Pettit test
#'
#' @param data_series A climatological time series.
#' @param ... Additional arguments to be passed into methods.
#'
#' @return Pettit test results.
#' @author Shadrack Kibet
#' @export
#'
#' @examples
pettit_test <- function(data_series, ...) {
  if(missing(data_series)) stop("argument 'data series' is missing with no default")
  n <- length(data_series)
  s2 <- 1:n
  s <- sapply(s2, function(k) {
    n <- length(data_series)
    r <- rank(data_series)
    b <- 0
    for (kk in 1:k) b <- b + r[kk]
    return(2 * b - kk * (n + 1))
  })
  s
}

#' Von Neumann's ratio test
#'
#' @param na.rm A logical value indicating whether missing values should be stripped before the computation proceeds.
#' @param data_series A climatological time series.
#' @param ... Additional arguments to be passed into methods.
#'
#' @return Von Neumann's ratio test results.
#' @author Shadrack Kibet
#' @export
#'
#' @examples
von_neumann_ratio_test <- function(data_series, na.rm = FALSE, ...) {
  if(missing(data_series)) stop("argument 'data series' is missing with no default")
  n <- length(data_series)
  m <- summary_mean(data_series, na.rm = na.rm)
  t1 <- data_series[1:(n - 1)]
  t2 <- data_series[2:n]
  if (summary_sum((data_series - m)^2, na.rm = na.rm) == 0) return(9999)
  N <- (summary_sum((t1 - t2)^2, na.rm = na.rm)) / summary_sum((data_series - m)^2, na.rm = na.rm)
  if (is.na(N)) return(9999)
  N
}

#' Bayesian test
#'
#' @param na.rm A logical value indicating whether missing values should be stripped before the computation proceeds.
#' @param data_series A climatological time series.
#' @param ... Additional arguments to be passed into methods.
#'
#' @return The Bayesian test results.
#' @author Shadrack Kibet
#' @export
#'
#' @examples
bayesian_test <- function(data_series, na.rm = FALSE, ...) {
  if(missing(data_series)) stop("argument 'data series' is missing with no default")
  dx <- summary_sd(data_series, na.rm = na.rm)
  n <- length(data_series)
  k <- 1:n
  sk <- sapply(k, function(kk) {
    m <- summary_mean(data_series, na.rm = na.rm)
    b <- summary_sum(data_series[1:kk] - m, na.rm = na.rm)
    return(b)
  })
  zk <- ((k * (n - k)^-0.5) * sk) / dx
  zk <- zk[!(is.nan(zk) | is.infinite(zk))]
  (A <- summary_sum(zk**2))
}

