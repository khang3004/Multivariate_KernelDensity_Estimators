source('local_linear_regression.R')
source(here::here('src/smoothing/local_linear_regression.R'))
#' Brute Force Leave-1-out Cross-validation to optimize H
#' @description
#' Calculate the CV score by physically removing 1 data point at a time,
#' retraining the model, and validating the removed point
#' This is the direct implementation
#' WARNING: Computationally expensive BigO(h*n^2)
#' @param x_train Numeric vector. Predictor
#' @param y_train Numeric vector. Response
#' @param h_grid Numeric vector. Grid of bandwidths to test.
#'
#' @return A list containing the optimal h and the CV score
#' @export
cv_brute_force <- function(x_train, y_train, h_grid,kernel = c("gauss", "epan")) {
  # 1. TYPE ASSERTION & VALIDATION
  stopifnot(length(x_train) == length(y_train))
  kernel <- match.arg(kernel)
  n <- length(x_train)
  # --- Worker function ---
  cv_brute_force_each_h <- function(h) {
    err <- 0 #Residual Sum of Squares
    # Loop through every data point (The Bottleneck)
    for (i in 1:n) {
      # 1. Blindfold: Remove the i-th point
      x_val <- x_train[i]   # The point held out
      y_val <- y_train[i]

      x_blind <- x_train[-i]  # The remaining n-1 points
      y_blind <- y_train[-i]

      # 2. Predict using the blindfolded model
      # We use our previously written local_linear_regression module
      # Note: This function is vectorized, but here we pass a single point x_val
      pred <- local_linear_regression(
        x_train = x_train[-i],
        y_train = y_train[-i],
        x_query = x_train[i],
        h = h,
        kernel = kernel
      )
      # 3. Accumulate Error
      # Handle NA (in case bandwidth is too small to find neighbors)
      if (!is.na(pred)) {
        err <- err + (y_val - pred)^2
      } else {
        next
      }
    }

    # Return Mean Squared Error (MSE)
    return(err / n)
  }

  #   3. Vectorization-based apply
  error_list <- vapply(h_grid, cv_brute_force_each_h, FUN.VALUE = numeric(1))

  #Find the index of h that minimize CV_error
  argmin_index <- which.min(error_list)

  # Return a LIST (Standard R practice for multiple return values)
  return(list(
    optimal_h = h_grid[argmin_index],
    min_cv_error = error_list[argmin_index],
    all_scores = error_list
  ))

}

#' 2. Shortcut LOOCV (Vectorized)
#' @description Implements Slide 91 formula using Vectorization.
#' Formula: Mean( ((y - y_hat) / (1 - w_ii))^2 )
#' @export

#' 2. Shortcut
#' @description Implements the following formula using Vectorization.
#' Formula: Mean( ((y - y_hat) / (1 - w_ii))^2 )
#' #' @param x_train Numeric vector. Predictor
#' @param y_train Numeric vector. Response
#' @param h_grid Numeric vector. Grid of bandwidths to test.
#'
#' @return A list containing the optimal h, minimum CV error,  the CV score
#' @export
cv_shortcut <- function(x_train, y_train, h_grid, kernel = c("gauss", "epan")) {
  stopifnot(length(x_train) == length(y_train))
  kernel <- match.arg(kernel)
  cv_shortcut_each_h <- function(h) {
    if (h <= 0) return(Inf)

    # 1. Fit Full Model (Get y_hat and W_ii)
    # Important: Pass x_query = x_train to get predictions at training points
    fit <- local_linear_regression(
      x_train = x_train,
      y_train = y_train,
      x_query = x_train,
      h = h,
      kernel = kernel,
      return_leverage = TRUE
    )

    y_hat <- fit$y_hat
    w_ii <- fit$w_ii

    # Safety Check
    if (any(is.na(y_hat))) return(Inf)

    # 2. Vectorized Calculation
    # Avoid division by zero if w_ii is close to 1 (Clipping to avoid NaN)
    w_ii[w_ii > 0.999] <- 0.999

    residuals <- y_train - y_hat
    penalized_residuals <- residuals / (1 - w_ii)

    # MSE
    return(mean(penalized_residuals^2))
  }

  error_list <- vapply(h_grid, cv_shortcut_each_h, FUN.VALUE = numeric(1))

  argmin_index <- which.min(error_list)
  return(list(optimal_h = h_grid[argmin_index], min_cv_error = error_list[argmin_index], all_scores = error_list))
}


#' 3. Generalized Cross-Validation (GCV)
#' @description
#' Implements GCV formula: Sum(Residuals^2) / (n - Sum(W_ii))^2
#' Optimized for speed.
#' @param x_train Numeric vector. Predictor
#' @param y_train Numeric vector. Response
#' @param h_grid Numeric vector. Grid of bandwidths to test.
#'
#' @return A list containing the optimal h, minimum CV error,  the CV score
#' @export
generalized_cv_shortcut <- function(x_train, y_train, h_grid, kernel = c("gauss", "epan")) {
  stopifnot(length(x_train) == length(y_train))
  kernel <- match.arg(kernel)
  n <- length(x_train)

  cv_gcv_each_h <- function(h) {
    if (h <= 0) return(Inf)

    fit <- local_linear_regression(
      x_train = x_train,
      y_train = y_train,
      x_query = x_train,
      h = h,
      kernel = kernel,
      return_leverage = TRUE
    )

    y_hat <- fit$y_hat
    w_ii <- fit$w_ii

    if (any(is.na(y_hat))) return(Inf)

    rss <- sum((y_train - y_hat)^2)

    trace_s <- sum(w_ii)

    # Denominator: (n - trace)^2
    denominator <- (n - trace_s)^2

    if (denominator < 1e-10) return(Inf)

    gcv_value <- rss / denominator
    return(gcv_value)
  }

  error_list <- vapply(h_grid, cv_gcv_each_h, FUN.VALUE = numeric(1))

  argmin_index <- which.min(error_list)
  return(list(optimal_h = h_grid[argmin_index],
              min_cv_error = error_list[argmin_index],
              all_scores = error_list))
}


#' 4. Grab CV \& Generalized CV into only function


#' Optimized & Robust Cross-Validation Selection
#'
#' @description
#' The ultimate bandwidth selector. Combines efficient computation (DRY)
#' with robust numerical handling for h -> 0 singularities.
#'
#' @param x_train Numeric vector.
#' @param y_train Numeric vector.
#' @param h_grid Numeric vector of bandwidths.
#' @param metric Character. "cv" (LOOCV), "gcv" (Generalized CV), or "both".
#' @param strategies Character vector. Strategies to handle singularity:
#'   - "bound": Enforce h > min_distance/2.
#'   - "clip":  Add epsilon to denominator (1 - w_ii) to prevent division by zero.
#'   - "robust": Use Median instead of Mean to aggregate LOOCV errors (ignore outliers).
#'
#' @return A list containing optimal h for selected metrics and a full score table.
#' @export
bw_cv_final <- function(x_train, y_train, h_grid,
                        kernel = c("gauss", "epan"),
                        metric = c("both", "cv", "gcv"),
                        strategies = c("bound", "clip")) {

  # --- 1. SETUP & VALIDATION ---
  stopifnot(length(x_train) == length(y_train))
  kernel <- match.arg(kernel)
  metric <- match.arg(metric)
  n <- length(x_train)

  # [STRATEGY 1: Bound]
  safe_h_limit <- 0
  if ("bound" %in% strategies) {
    min_dx <- min(diff(sort(unique(x_train))))
    safe_h_limit <- min_dx * 0.51
  }

  # --- 2. WORKER FUNCTION ---
  calc_single_h <- function(h) {
    # A. Pre-check Safety
    if (h < safe_h_limit) {
      return(c(cv = Inf, gcv = Inf, status = "Skipped (Too Small)"))
    }

    fit <- local_linear_regression(
      x_train, y_train, x_query = x_train,
      h = h,
      kernel = kernel, 
      return_leverage = TRUE
    )

    y_hat <- fit$y_hat
    w_ii  <- fit$w_ii

    # Check Core Failure
    if (any(is.na(y_hat))) {
      return(c(cv = Inf, gcv = Inf, status = "Failed (Core NaN)"))
    }

    residuals <- y_train - y_hat
    res_cv <- NA
    res_gcv <- NA
    status <- "OK"

    # C. Calculate LOOCV
    if (metric %in% c("cv", "both")) {
      denom_cv <- 1 - w_ii

      # [STRATEGY 2: Clip]
      if ("clip" %in% strategies) {
        denom_cv <- pmax(denom_cv, 1e-6) # Epsilon trick
      } else {
        if (any(abs(denom_cv) < 1e-10)) status <- "Singularity Warning"
      }

      sq_errors <- (residuals / denom_cv)^2

      # [STRATEGY 3: Robust Aggregation]
      if ("robust" %in% strategies) {
        res_cv <- median(sq_errors, na.rm = TRUE)
      } else {
        res_cv <- mean(sq_errors, na.rm = TRUE)
      }
    }

    # D. Calculate GCV
    if (metric %in% c("gcv", "both")) {
      trace_s <- sum(w_ii)
      denom_gcv <- (1 - trace_s / n)^2

      if (denom_gcv < 1e-10) denom_gcv <- 1e-10

      res_gcv <- sum(residuals^2) / (n * denom_gcv)
    }

    return(c(cv = res_cv, gcv = res_gcv, status = status))
  }

  # --- 3. VECTORIZED EXECUTION ---
  raw_results <- lapply(h_grid, calc_single_h)

  results_df <- do.call(rbind, lapply(raw_results, function(x) {
    data.frame(
      h = NA, # Placeholder
      cv = as.numeric(x["cv"]),
      gcv = as.numeric(x["gcv"]),
      status = x["status"],
      stringsAsFactors = FALSE
    )
  }))
  results_df$h <- h_grid

  # --- 4. FIND OPTIMAL H ---
  opt_h_cv <- NA
  opt_h_gcv <- NA

  # Helper find minimum & skip NAs
  find_best <- function(scores) {
    valid <- is.finite(scores)
    if (!any(valid)) return(NA)
    idx <- which(scores == min(scores[valid]))[1]
    return(h_grid[idx])
  }

  if (metric %in% c("cv", "both")) opt_h_cv <- find_best(results_df$cv)
  if (metric %in% c("gcv", "both")) opt_h_gcv <- find_best(results_df$gcv)

  return(list(
    opt_h_cv = opt_h_cv,
    opt_h_gcv = opt_h_gcv,
    scores = results_df,
    metric_used = metric,
    strategies_used = strategies
  ))
}
