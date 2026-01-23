library(testthat)

# 1. Source the modules (Adjust path relative to this file)
# Assuming running from project root or having set wd correctly
library(here)
source(here::here("src", "smoothing", "local_linear_regression.R"))
source(here::here("src", "smoothing", "bandwidth_smoothing_selection", "cross_validation.R"))

# --- TEST GROUP 1: BASIC FUNCTIONALITY ---
test_that("bw_cv_final returns correct structure", {
  # Mock data: Simple linear relation y = 2x
  x <- 1:10
  y <- 2 * x
  h_grid <- c(0.5, 1.0, 2.0)

  # Run function
  res <- bw_cv_final(x, y, h_grid, metric = "both")

  # Assertions
  expect_type(res, "list")
  expect_named(res, c("opt_h_cv", "opt_h_gcv", "scores", "metric_used", "strategies_used"))
  expect_s3_class(res$scores, "data.frame")
  expect_equal(nrow(res$scores), length(h_grid))
})

# --- TEST GROUP 2: EDGE CASES & SINGULARITY ---
test_that("bw_cv_final handles h -> 0 correctly (Singularity)", {
  # Scenario: Two points extremely close -> causing W_ii ~ 1 if h is small
  x_nasty <- c(1.0, 1.00001, 2.0, 3.0)
  y_nasty <- c(1, 1, 2, 3)

  # A bandwidth smaller than the gap (dangerous!)
  h_dangerous <- 0.000005
  h_safe <- 0.5
  h_grid <- c(h_dangerous, h_safe)

  # Case A: Strategy "bound" should SKIP the dangerous h
  res_bound <- bw_cv_final(x_nasty, y_nasty, h_grid, strategies = c("bound"))
  status_dangerous <- res_bound$scores[res_bound$scores$h == h_dangerous, "status"]
  expect_match(status_dangerous, "Skipped", info = "Should skip h < safe_limit")

  # Case B: Strategy "clip" should CALCULATE but return FINITE score (not Inf)
  # We disable 'bound' to force it to calculate
  res_clip <- bw_cv_final(x_nasty, y_nasty, h_grid, strategies = c("clip"))
  cv_score <- res_clip$scores[res_clip$scores$h == h_dangerous, "cv"]

  expect_true(is.finite(cv_score), info = "Clipping should prevent Inf")
  expect_gt(cv_score, 1000, info = "Score should be very high due to penalty")
})

# # --- TEST GROUP 3: RIDGE REGULARIZATION ---
# test_that("Ridge regularization works without crashing", {
#   x <- runif(20)
#   y <- sin(x) + rnorm(20)
#   h_grid <- c(0.1)
#
#   # Run with Ridge strategy
#   res <- bw_cv_final(x, y, h_grid, strategies = c("ridge"), lambda = 1e-4)
#
#   expect_equal(res$scores$status[1], "OK")
#   expect_false(any(is.na(res$scores$cv)))
# })

# --- TEST GROUP 4: GCV VS CV ---
test_that("GCV and CV produce different but correlated results", {
  x <- seq(0, 10, length.out = 50)
  y <- sin(x) + rnorm(50, sd = 0.2)
  h_grid <- seq(0.5, 2.0, by = 0.5)

  res <- bw_cv_final(x, y, h_grid, metric = "both")

  # Ensure both are computed
  expect_false(all(is.na(res$scores$cv)))
  expect_false(all(is.na(res$scores$gcv)))

  # They shouldn't be exactly equal (mathematically impossible usually)
  expect_false(all(res$scores$cv == res$scores$gcv))
})