
# Homogenization tests ----------------------------------------------------

#' Buishand range test
#'
#' @param data
#' @param series
#' @param na.rm
#'
#' @return
#' @export
#'
#' @examples
Buishand_range_test <- function(data, series = "", na.rm = FALSE) {
  if(missing(series)) stop("Data series name is missing")
  data_series <- data[,series]
  n <- length(data_series)
  s2 <- 1:n
  s <- sapply(s2, function(k) {
    m <- summary_mean(data_series, na.rm = na.rm)
    b <- summary_sum(data_series[1:k] - m, na.rm = na.rm)
    return(b)
  })
  s / summary_sd(data_series, na.rm = na.rm) / sqrt(n)
}

#' craddock test
#'
#' @param data
#' @param col_name1
#' @param col_name2
#' @param na.rm
#' @param dp
#' @param ...
#'
#' @examples
Craddock_test <- function(data, col_name1 = "", col_name2 = "", na.rm = FALSE, dp = 2, ...) {
  if (missing(col_name1)||missing(col_name2)) stop("Both columns names must be supplied")
  data <- as.data.frame(data)
  col1 <- data[, col_name1]
  col2 <- data[, col_name2]
  mean_col1 <- summary_mean(col1, na.rm = na.rm)
  mean_col2 <- summary_mean(col2, na.rm = na.rm)
  c2C <- col2 + mean_col1 - mean_col2
  col_diff <- col1 - c2C
  round(cumsum(col_diff), digits = dp)
}

#' Standard Normal Homogeneity test(SNHT)
#'
#' @param data
#' @param series
#' @param na.rm
#' @param dp
#' @param ...
#'
#' @examples
snht <- function(data, series = "", na.rm = FALSE, dp = 2, ...) {
  if(missing(series)) stop("Data series name is missing")
  data <- as.data.frame(data)
  data_series <- data[,series]
  tsn <- (data_series - summary_mean(data_series, na.rm = na.rm)) / summary_sd(data_series, na.rm = na.rm)
  n <- length(data_series)
  zvec <- vector("numeric", n - 1)
  for (v in 1:(n - 1)) {
    zvec[v] <- v * (summary_mean(tsn[1:v], na.rm = na.rm))^2 +
      (n - v) * (summary_mean(tsn[(v + 1):n], na.rm = na.rm))^2
  }
  round(zvec, digits = dp)
}

#' pettit test
#'
#' @param data
#' @param series
#' @param ...
#'
#' @examples
Pettitt_test <- function(data, series = "", ...) {
  if(missing(series)) stop("Data series name is missing")
  data <- as.data.frame(data)
  data_series <- data[,series]
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

#' Von Neumman ratio test
#'
#' @param data
#' @param series
#' @param na.rm
#' @param ...
#'
#' @examples
von_neumann_ratio_test <- function(data, series = "", na.rm = FALSE, ...) {
  if(missing(series)) stop("Data series name is missing")
  data_series <- data[,series]
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
#' @param data
#' @param series
#' @param na.rm
#' @param ...
#'
#' @examples
bayesian_test <- function(data, series = "", na.rm = FALSE, ...) {
  if(missing(series)) stop("Data series name is missing")
  data_series <- data[,series]
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

