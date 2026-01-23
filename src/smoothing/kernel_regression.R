source("src/kernel_density_estimation/gaussian_kernel.R")
#' Nadaraya - WATSON kernel Regression
#'
#' @description
#' Implements the Nadaraya-Watson estimator (Local Constant Smoothing).
#' Formula: y_hat = Sum (W_i*Y_iz)/sum(W_i)
#'
#' @param x_train NUmeric vector. Observed predictor values.
#' @param y_train Numeric vector. Observed reponse values.
#' @param x_query Numeric vector. Points to estimate the trend.
#' @param h numeric value. Bandwidth (smoothing parameter).
#' @return Numeric vector of estimated values.
#' @export
nadaraya_watson <- function(x_train, y_train, x_query, h = 0.3, kernel_type = "gaussian") {

  # ---1. TYPE ASSERTION SAFETY & VALIDATION ---
  stopifnot(is.numeric(x_train), is.numeric(y_train))
  stopifnot(length(x_train) == length(y_train))  #check data length matching
  stopifnot(is.numeric(h), h > 0)   # Check if Bandwidth is positive?

  # ---2. WORKER FUNCTION (Single Point Estimation)
  estimate_single <- function(x_q) {
    kernel_weights <- compute_contribution_gaussian(x = x_q, observation = x_train, h = h)
    numerator <- sum(y_train * kernel_weights)
    wi_sum <- sum(kernel_weights)
    return(numerator / wi_sum)
  }

  # ---3. VECTORIZATION (Performance Optimization)
  # Apply the estimator to all query points
  #vapply guarantees the ouput is a numeric vector
  y_hat <- vapply(x_query, estimate_single, FUN.VALUE = numeric(1))

  return(y_hat)
}

# vectoried_nadaraya_watson <- Vectorize(FUN = nadaraya_watson, vectorize.args = "x_query")






# _____________________________________-
#' Funstion & Modules in the lECTURE 2
kernel_gauss <- function(u) exp(-u^2 / 2) / sqrt(2 * pi)

kernel_reg <- function(x, y, x_eval, h, kernel = "gauss") {
  wi <- kernel_gauss((x_eval - x) / h)
  wi_sum <- sum(wi)
  w <- wi / wi_sum
  y_hat <- sum(y * w)
  return(y_hat)
}

kernel_reg <- Vectorize(FUN = kernel_reg, vectorize.args = "x_eval")


# CV_1 ----
cv1 <- function(x, y, h, kernel = "gauss") {
  n <- length(x)
  err <- numeric(n)
  for (i in 1:n) {
    x_new <- x[-i]
    y_new <- y[-i]
    y_i <- kernel_reg(
      x = x_new,
      y = y_new,
      x_eval = x[i],
      h = h,
      kernel = kernel
    )
    err[i] <- (y[i] - y_i)^2
  }

  cv_hat <- mean(err)
  return(cv_hat)
}

##CV----

cv1_h <- Vectorize(FUN = cv1, vectorize.args = "h")


cv2_h <- function(x, y, h, kernel = "gauss") {
  n <- length(x)
  w <- matrix(0, n, n)
  wi0 <- numeric(n)
  for (i in 1:n) {
    wi0[i] <- sum(kernel_gauss((x[i] - x) / h))
  }
  w0 <- kernel_gauss(0) / wi0

  m0 <- w0 * sum(y)
  return(m0)
}
  