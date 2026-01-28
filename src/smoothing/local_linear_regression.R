#' Kernel function for local polynomial regression
#'
#' Computes kernel weights used in local linear regression.
#' The Gaussian kernel is used in unnormalized form, since any
#' multiplicative constant cancels out in the estimator.
#'
#' @param u Numeric vector of scaled distances, \eqn{u = (x - x_0)/h}.
#' @param kernel Character string specifying the kernel:
#'   \code{"gauss"} (Gaussian, unnormalized) or \code{"epan"} (Epanechnikov).
#'
#' @return Numeric vector of kernel weights.
#'
#' @export

kernel_fun <- function(u, kernel = c("gauss", "epan")) {
  kernel <- match.arg(kernel)
  
  if (kernel == "gauss") {
    # Gaussian kernel (unnormalized): K(u) = exp(-u^2 / 2)
    exp(-0.5 * u^2)
  } else {
    # Epanechnikov kernel: K(u) = 3/4 (1 - u^2) I(|u| <= 1)
    0.75 * (1 - u^2) * (abs(u) <= 1)
  }
}

#' Local Linear Estimator (Enhanced with Diagnostics)
#'
#' @description
#' Implements Local Linear Regression.
#' Updated to support returning Leverage (Hat Matrix Diagonal) for CV Shortcut.
#'
#' @param x_train Numeric vector. Observed data X.
#' @param y_train Numeric vector. Observed data Y.
#' @param x_query Numeric vector. Points to estimate.
#' @param h Numeric value. Bandwidth.
#' @param kernel Character. Kernel type ("gauss" or "epan").
#' @param return_leverage Logical. If TRUE, returns a list containing y_hat and W_ii.
#'                        Only valid when x_query is practically same as x_train (for CV).
#'
#' @export
local_linear_regression <- function(
  x_train, y_train, x_query, h,
  kernel = c("gauss", "epan"),
  return_leverage = FALSE
) {
  # --- 1. VALIDATION ---
  stopifnot(length(x_train) == length(y_train), h > 0)
  kernel <- match.arg(kernel)
  # Kernel value at zero, consistent with kernel_fun
  # Gaussian: K(0) = 1 (unnormalized)
  # Epanechnikov: K(0) = 0.75
  K_0 <- if (kernel == "gauss") 1 else 0.75
  

  # --- 2. WORKER FUNCTION (Returns LIST to handle multiple outputs) ---
  estimate_single_detailed <- function(x0) {
    # a. Calculate Diff & Weights
    diffs <- x_train - x0
    weights <- kernel_fun(diffs / h, kernel)

    # Threshold checks
    if (sum(weights) < 1e-10) return(list(pred = NA, leverage = NA))

    # b. Moments S0, S1, S2 (Omitted 1/n)
    s0 <- sum(weights)
    s1 <- sum(weights * diffs)
    s2 <- sum(weights * diffs^2)

    # c. Prediction Logic (y_hat)
    sum_wy <- sum(weights * y_train)
    sum_w_diff_y <- sum(weights * diffs * y_train)

    det <- s2 * s0 - s1^2

    if (abs(det) < 1e-10) return(list(pred = NA, leverage = NA))

    y_hat <- (s2 * sum_wy - s1 * sum_w_diff_y) / det

    # d. Leverage Logic (W_ii) - Only calculated if requested
    # Formula: W_ii = (K(0) * S2) / Determinant
    # Note: This technically assumes x0 is one of the training points (distance -> 0)
    # But mathematically it represents the "weight at the center" regardless.
    leverage <- NA
    if (return_leverage) {
      leverage <- (K_0 * s2) / det
      # Cap leverage to avoid singular LOO error
      if (leverage > 0.999) leverage <- 0.999
    }

    return(list(pred = y_hat, leverage = leverage))
  }

  # --- 3. VECTORIZATION ---
  # Since we might return a complex object, standard vapply is tricky.
  # We use lapply and then bind results.

  results_list <- lapply(x_query, estimate_single_detailed)

  # Extract y_hat vector
  y_hat <- vapply(results_list, function(r) r$pred, numeric(1))

  # --- 4. RETURN logic ---
  if (return_leverage) {
    w_ii <- vapply(results_list, function(r) r$leverage, numeric(1))
    return(list(y_hat = y_hat, w_ii = w_ii))
  } else {
    return(y_hat)
  }
}