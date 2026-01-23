source('./src/smoothing/local_polynomial_regression.R')
marsbig <- read.table('./data/marsbig.dat', header = TRUE)

#' Local Polynomial Regression Examples 2
#' @param Planet Radius
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
  x_train <- marsbig_o$radius
  y_train <- marsbig_o$temperature
  x_query <- seq(min(x_train), max(x_train), length.out = 500)
  
  # Use local polynomial regression with p = 1, 2; h = NULL
  y_hat_1 <- local_polynomial_regression(x_train, y_train, x_query, p = 1)
  y_hat_2 <- local_polynomial_regression(x_train, y_train, x_query, p = 2)
  
  # Plot sample training and sample predicting
  plot(
    x_train, y_train,
    pch = 16, col = "gray",
    xlab = "Planet Radius",
    ylab = "Temperatures",
    main = paste0("Orbit ", o)
  )
  
  lines(x_query, y_hat_1, col = 'blue', lwd = 2)
  lines(x_query, y_hat_2, col = 'red', lwd = 2)
}
# Summary: The temperature decreases as the radius increases from 3375 to 3400, and then begins
# to increase gradually at orbit 1, 2, 6. The temperature increases slightly and then 
# gradually decreases as the radius increases at orbit 3, 4, 5, 7.
