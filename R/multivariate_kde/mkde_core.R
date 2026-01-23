source("R/multivariate_kde/mkde_kernels.R")
source("R/multivariate_kde/bw_selection_multivariate.R")
#' @description
#' Computes density estimate using Product Kernel formulation.
#' High optimized using matrix operations instead of loops.
#'
#' Assume all dimension is inpendent
#' @parameter data Numeric matrix (n \times p). Training data.
#' @parameter x_eval Numeric matrix (m\times p) or vector (p). Points to evaluate density
#' @parameter h Numeric vector (p). Bandwidth vector.
#' @return numeric vector of densities
#' @export

product_kernel_mkde <- function(x_eval, data, h) {
  # --- 1. PRE-PRECESSING ---
  # --- TYPE SAFETY ---
  if (is.data.frame(data)) data <- as.matrix(data)
  stopifnot(is.matrix(data))
  # Handle x_eval being a single vector -> convert to 1-row matrix
  if (is.vector(x_eval)) {
    x_eval <- matrix(x_eval, nrow = 1)
  }
  if (is.data.frame(x_eval)) x_eval <- as.matrix(x_eval)

  # Get dimensions of data matrix
  n <- nrow(data)
  p <- ncol(data)
  stopifnot(ncol(x_eval) == p) #Ensure dimensions match

  # Auto-select bandwidth if missing
  if (is.null(h)) h <- multi_silverman_bandwidth(data)
  stopifnot(length(h) == p)

  # ---2. Estimation (Worker function)----
  # Calculate density for a SINGLE query point (x_q) using Vectorization
  estimate_each_point <- function(x_q) {
    # x_q is a vector of length p

    # A. Broadcasting Difference:
    # Create matrix of differences (n\times p)
    # We use t(t(data)-x_q) to subtract vector x_q form every row of data
    diffs <- t(t(data) - x_q)

    #B. Scale by h (Broadcasting division)
    scaled_diffs <- t(t(diffs) / h)

    # C. Apply Univariate Kernel to All elements
    # Result is (n\times p) matrix of kernel values
    k_vals <- gaussian_kernel(scaled_diffs)

    # D. Product across dimensions (The "Product in Product Kernel)
    # Result is vector of length n (Weight of each observation)
    w_i <- apply(k_vals, 1, prod)

    # E. Sum & Normaliza
    # Density = (1/(n*prod(h)))*Sum(w_i)
    density <- sum(w_i) / (n * prod(h))
    return(density)
  }


  # --- 3. Execution ---
  # Apply to every row in x_eval
  # 'apply' iterates over rows (MARGINS=1)
  y_hat <- apply(x_eval,1,estimate_each_point) # a matrix has ncols = p
  return(y_hat)
}
