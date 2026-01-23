#' Scott's Rule of Thumb for Multivariate KDE Bandwidth
#' @description
#' Calculate the optimal bandwidth for each dimension using Scott's Rule
#' Assumes the underlying density p(x) is somewhat Gaussian-like
#' Formula: h_j = (4/(n*(p+2)))^{1/(p+4)}  \sigma_hat_i
#' @parameter X_train: (n*p)-dim matrix
#' @return h: a p-dim vector
#' @export
#'
multi_silverman_bandwidth <- function(data) {
  # --- TYPE SAFETY ---
  if (is.data.frame(data)) data <- as.matrix(data)
  stopifnot(is.matrix(data))

  # Get dimensions of data matrix
  n <- nrow(data)
  p <- ncol(data)

  #Calculate Standard Deviation for each dimension
  # apply(data, 2, sd) means appply sd function to columns (dim 2)


  sigma_vector <- apply(data, 2, sd)
  bandwidth_vector <- (4 / (n * (p + 2)))^{ 1 / (p + 4) } * sigma_vector
  return(bandwidth_vector)
}
