#' Data Whitening (Sphering) Module
#'
#' @description
#' Transforms correlated data into uncorrelated (white) data using Mahalanobis transformation.
#' Essential for MKDE when using Product Kernels on correlated datasets.

#' Compute Whitening Matrix
#' @param x Numeric matrix (n x d).
#' @return A list containing:
#'   - mean: Vector of means.
#'   - W: The whitening matrix (d x d). \Sigma^(-1/2)
#'   - W_inv: The inverse matrix (d x d) for mapping back. \Sigma^(1/2)
get_whitening_params <- function(x) {
  # 1. Calculate Statistics
  mu <- colMeans(x)
  Sigma <- cov(x)

  # 2. Eigen Decomposition (Robust way to compute matrix power)
  # Sigma = U * D * U'
  eigen_decomp <- eigen(Sigma)

  U <- eigen_decomp$vectors
  D <- eigen_decomp$values

  # Handle numerical instability (negative eigenvalues close to 0)
  D <- pmax(D, 1e-10)

  # 3. Compute W = Sigma^(-1/2) = U * D^(-1/2) * U'
  # D^(-1/2) is just 1/sqrt(eigenvalues) on diagonal
  D_inv_sqrt <- diag(1 / sqrt(D))

  W <- U %*% D_inv_sqrt %*% t(U)

  # Compute Inverse (for reference/plotting) if needed
  # W_inv = Sigma^(1/2)
  D_sqrt <- diag(sqrt(D))
  W_inv <- U %*% D_sqrt %*% t(U)

  return(list(mean = mu, W = W, W_inv = W_inv))
}

#' Apply Whitening Transformation
#' @param x Numeric matrix.
#' @param params List from get_whitening_params.
#' @return Transformed matrix Z where cov(Z) approx Identity.
whiten_data <- function(x, params) {
  # Center the data: X - mu
  x_centered <- sweep(x, 2, params$mean, "-")

  # Rotate and Scale: (X - mu) %*% W
  # Note: Standard algebra is W * x_vec, but for data matrix NxD, it's X %*% W
  z <- x_centered %*% params$W

  return(z)
}