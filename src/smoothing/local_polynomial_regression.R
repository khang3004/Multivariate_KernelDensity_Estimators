source(here::here('src/smoothing/bandwidth_smoothing_selection/cross_validation.R'))

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
  stopifnot(length(x_train) == length(y_train), p >= 0)
  
  # Check if bandwidth is Null
  if (is.null(h)) {
    # Define grid of bandwidth
    # min_dist: min distance between each pair of x_train
    # range_x: Diff between max and min of x_train
    # Note: range_x / 2 avoid over-smoothing
    min_dist <- min(diff(sort(unique(x_train))))
    range_x <- diff(range(x_train))
    grid <- seq(min_dist, range_x / 2, length.out = 20)
    
    # If length of x_train <= 10 ^ 3: using metric = 'cv', else: using metric = 'gcv'
    n <- length(x_train)
    if (n <= 10 ^ 3) {
      fit <- bw_cv_final(x_train, y_train, grid, metric = 'cv')
      h <- fit$opt_h_cv
    } else {
      fit <- bw_cv_final(x_train, y_train, grid, metric = 'gcv')
      h <- fit$opt_h_gcv
    }
  } else {
    if (h <= 0) {
      stop('h must be positive or null')
    }
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
    
    # Implement Ridge Regularization to handle singular matrix
    lambda <- 10 ^ (-6)
    I <- diag(p + 1)
    
    # Calculate y_hat by formula e^T * (Z^T * W * Z + lambda * I)^(-1) * W * Y; 
    # e = (1, 0, 0,.., 0)
    e <- c(1, rep(0, p))
    y_hat = t(e) %*% solve(t(Z) %*% W %*% Z + lambda * I) %*% t(Z) %*% W %*% y_train
    
    return(y_hat = y_hat)
  }
  
  # --- 3. VECTORIZATION ---
  y_hat <- vapply(x_query, estimate_single_detailed, FUN.VALUE = numeric(1))
  
  # --- 4. RETURN logic ---
  return(y_hat = y_hat)
}
