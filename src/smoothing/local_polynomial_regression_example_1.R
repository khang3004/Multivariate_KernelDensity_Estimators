source(here::here('src/smoothing/local_polynomial_regression.R'))
mars <- read.table(here::here("data", "mars.dat"), header = TRUE)

#' Local Polynomial Regression Examples 1
#' @description Predict temperature trend based on the changes of planet radius.
#' 
#' @param Planet Radius
#' @param Temperatures
#' 
#' @return Temperature trend with variable radius

# Get sample training, sample predicting
x_train <- mars$radius
y_train <- mars$temperature
x_query <- seq(min(x_train), max(x_train), length.out = 200)

# Use local polynomial regression with p = 1, 2; h = NULL
y_hat_1 <- local_polynomial_regression(x_train, y_train, x_query, p = 1)
y_hat_2 <- local_polynomial_regression(x_train, y_train, x_query, p = 2)

# Plot sample training and sample predicting
plot(
  x_train, y_train,
  pch = 16, col = "gray",
  xlab = "Planet Radius",
  ylab = "Temperatures",
  main = "Temperature trend with variable radius"
)

lines(x_query, y_hat_1, col = 'blue', lwd = 2)
lines(x_query, y_hat_2, col = 'red', lwd = 2)

legend("topright", legend=c("p=1", "p=2"), col=c("blue", "red"), lty=1)

