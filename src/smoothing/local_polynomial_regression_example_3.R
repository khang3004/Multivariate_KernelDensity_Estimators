source('./src/smoothing/local_polynomial_regression.R')
marsbig <- read.table('./data/marsbig.dat', header = TRUE)

marsbig

#' Local Polynomial Regression Examples 3
#' @param Pressure
#' @param Temperatures
#' @param Orbit
#' @return Temperature trend with variable radius each orbit

# Get unique orbit
orbits <- unique(marsbig$orbit)

par(mfrow = c(2, 4))
for (o in orbits) {
  # Extract sample based on orbit
  marsbig_o <- subset(marsbig, orbit == o)
  
  # Get sample training, sample predicting
  x_train <- marsbig_o$pressure
  y_train <- marsbig_o$temperature
  x_query <- seq(min(x_train), max(x_train), length.out = 200)
  
  # Use local polynomial regression with p = 1, 2; h = NULL
  y_hat_1 <- local_polynomial_regression(x_train, y_train, x_query, p = 1)
  y_hat_2 <- local_polynomial_regression(x_train, y_train, x_query, p = 2)
  
  # Plot sample training and sample predicting
  plot(
    x_train, y_train,
    pch = 16, col = "gray",
    xlab = "Pressure",
    ylab = "Temperatures",
    main = paste0("Orbit ", o)
  )
  
  lines(x_query, y_hat_1, col = 'blue', lwd = 2)
  lines(x_query, y_hat_2, col = 'red', lwd = 2)
}
# Summary: Error system is computationally singular: reciprocal condition number = 2.95886e-17
# This means that the matrix is nearly singular or these two feature have a strong correlation
