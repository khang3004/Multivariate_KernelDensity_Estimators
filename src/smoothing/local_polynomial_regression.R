source('./src/smoothing/bandwidth_smoothing_selection/cross_validation.R')

#' Local Polynomial Regression
#' @description Implements the Local Polynomial Regression
#' 
#' @param x_train: Numeric vector of observed sample X
#' @param y_train: Numeric vector of observed sample Y
#' @param x_query: Numeric vector of points
#' @param p: integer (p = 0, 1, 2,...)
#' @param h: bandwidth (default = Null)
#' 
#' @return y_hat: Numeric vector of predicted sample Y

local_polynomial_regression <- function(x_train, y_train, x_query, p, h = NULL) {
  # --- 1. VALIDATION ---
  stopifnot(length(x_train) == length(y_train), h > 0, p >= 0)
  
  # Check bandwidth
  if (is.null(h)) {
    fit <- generalized_cv_shortcut(x_train, y_train, seq(2, 100, 2))
    h <- fit$optimal_h
  }
  # Internal Kernel: K(u) = 1 / sqrt(2 * pi) * exp(-0.5 * u^2)
  .gaussian_kernel <- function(u) 1 / sqrt(2 * pi) * exp(-0.5 * u ^ 2)
  
  # --- 2. WORKER FUNCTION (Returns LIST to handle multiple outputs) ---
  estimate_single_detailed <- function(x0) {
    # a. Calculate Diff & Weights
    diffs <- x_train - x0
    weights <- .gaussian_kernel(diffs / h)
    
    # Threshold checks
    if (sum(weights) < 1e-10) return(pred = NA)
    
    # Create weight matrix W
    W <- diag(weights)
    
    # Create design matrix Z
    Z <- sapply(0:p, function(u) diffs ^ u)
    
    # Calculate y_hat by formula e^T * (Z^T * W * Z)^(-1) * W * Y; e = (1, 0, 0,.., 0)
    e <- c(1, rep(0, p))
    y_hat = t(e) %*% solve(t(Z) %*% W %*% Z) %*% t(Z) %*% W %*% y_train
    
    return(y_hat = y_hat)
  }
  
  # --- 3. VECTORIZATION ---
  y_hat <- lapply(x_query, estimate_single_detailed)
  
  # --- 4. RETURN logic ---
  return(y_hat = y_hat)
}
