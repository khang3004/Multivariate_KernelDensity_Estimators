set.seed(1)
#5
#a

# ham m(x) = 1 + sin(x^2)/x^2
m <- function(x){
  1 + sin((x^2))/x^2
}

# X pp deu tren (0,pi) : Uniform(0,pi)
X <- runif(100,0,pi)

# Độ lệch chuẩn của ei: sqrt(m(x)^2 / 64)
sd_eps <- sqrt(m(X)^2 / 64)

eps <- rnorm(100,0,sd_eps)

Y <-m(X) + eps

#b
kernel_fun <- function(u, kernel){
  if (kernel == "gauss") {
    return(exp(-u^2/2)/sqrt(2*pi))
  }
  else if (kernel == "epan") {
    return(0.75*(1 - u^2)*(abs(u) <= 1))
  }
  else stop("Unknown kernel")
}

#' Local Linear Regression
#'
#' @description
#' Performs local linear regression (kernel smoothing) at given evaluation points.
#' Formula: y_hat = mean(ww * y), where ww are local linear weights.
#'
#' @param x Numeric vector. Observed predictor values.
#' @param y Numeric vector. Observed response values.
#' @param x_eval Numeric vector. Points at which to evaluate the local linear regression.
#' @param h Numeric. Bandwidth (smoothing parameter).
#' @param kernel Character. Kernel type: "gauss" for Gaussian, "epan" for Epanechnikov.
#'
#' @return Numeric vector of estimated values at `x_eval`.
#' @export
loclin_reg <- function(x, y, x_eval, h, kernel){
  u <- x - x_eval
  uh <- u/h
  w <- kernel_fun(uh,kernel)
  s0 <- mean(w)
  s1 <- mean(w * (x - x_eval))
  s2 <- mean(w * (x - x_eval)^2)
  ww <- ((s2 - s1 * (x - x_eval)) * w)/(s2 * s0 - s1^2)
  y_hat <- mean(ww * y)
  return(y_hat)
}

loclin_reg <- Vectorize(FUN = loclin_reg, vectorize.args = "x_eval")

#' Rule-of-Thumb Bandwidth for Local Linear Regression
#'
#' @description
#' Computes the optimal bandwidth using a plug-in Rule-of-Thumb approach,
#' based on a quartic polynomial fit of the data.
#'
#' @param x Numeric vector. Observed predictor values.
#' @param y Numeric vector. Observed response values.
#' @param kernel Character. Kernel type: "gauss" or "epan".
#'
#' @return Numeric. Estimated bandwidth h_RT.
#' @export
h_ll_RT <- function(x, y, kernel) {
  if (kernel == "gauss") {
    R_k <- 0.5/sqrt(pi)
    sigma2_k <- 1
  } 
  else if (kernel == "epan") {
    R_k <- 3/5
    sigma2_k <- 1/5
  } 
  else stop("Unknown kernel")
  # Quartic fit
  mod_Q <- lm(y ~ poly(x, raw = TRUE, degree = 4))
  # Estimates of unknown quantities
  int_sigma2_hat <- diff(range(x))*sum(mod_Q$residuals^2)/mod_Q$df.residual
  m_2prime <- (2 * mod_Q$coefficients[3] + 6 * mod_Q$coefficients[4] * x +
                 12 * mod_Q$coefficients[5] * x^2)^2
  a_hat <- mean(m_2prime)
  # h_RT
  h_RT <- ((R_k*int_sigma2_hat)/(a_hat*sigma2_k*length(x)))^(1/5)
  return(h_RT)
}

#' Cross-Validation Bandwidth for Local Linear Regression
#'
#' @description
#' Computes the optimal bandwidth by leave-one-out cross-validation (LOOCV).
#'
#' @param x Numeric vector. Observed predictor values.
#' @param y Numeric vector. Observed response values.
#' @param h Numeric vector. Candidate bandwidth values.
#' @param kernel Character. Kernel type: "gauss" or "epan".
#'
#' @return Numeric vector. CV score for each bandwidth in `h`.
#' @export
CV_loclin <- function(x, y, h, kernel) {
  n <- length(y)
  S_diag <- numeric(n)
  for (i in 1:n) {
    u <- x - x[i]
    w <- kernel_fun(u/h, kernel)
    a0 <- mean(w)
    a1 <- mean(w * u)
    a2 <- mean(w * u^2)
    S_diag[i] <- (a2 * w[i])/(a2 * a0 - a1^2)/n
  }
  y_hat <- loclin_reg(x = x, y = y, x_eval = x, h = h, kernel = kernel)
  loocv <- mean(((y - y_hat)/(1 - S_diag))^2)
  return(loocv)
}

CV_loclin <- Vectorize(FUN = CV_loclin, vectorize.args = "h")

#' Generalized Cross-Validation Bandwidth for Local Linear Regression
#'
#' @description
#' Computes the optimal bandwidth using generalized cross-validation (GCV).
#'
#' @param x Numeric vector. Observed predictor values.
#' @param y Numeric vector. Observed response values.
#' @param h Numeric vector. Candidate bandwidth values.
#' @param kernel Character. Kernel type: "gauss" or "epan".
#'
#' @return Numeric vector. GCV score for each bandwidth in `h`.
#' @export
GCV_loclin <- function(x, y, h, kernel) {
  n <- length(y)
  S_diag <- numeric(n)
  for (i in 1:n) {
    u <- x - x[i]
    w <- kernel_fun(u/h, kernel)
    a0 <- mean(w)
    a1 <- mean(w * u)
    a2 <- mean(w * u^2)
    S_diag[i] <- (a2 * w[i])/(a2 * a0 - a1^2)/n
  }
  y_hat <- loclin_reg(x = x, y = y, x_eval = x, h = h, kernel = kernel)
  gcv <- sum((y - y_hat)^2)/(n - sum(S_diag))^2
  return(gcv)
}

GCV_loclin <- Vectorize(FUN = GCV_loclin, vectorize.args = "h")

h_grid <- seq(0.05, 1, length.out = 100)
kernels <- c("gauss", "epan")

res_list <- list()

for (ker in kernels) {
  # Rule of Thumb
  h_rt <- h_ll_RT(X, Y, kernel = ker)
  
  # CV
  cv_vals <- CV_loclin(X, Y, h = h_grid, kernel = ker)
  h_cv <- h_grid[which.min(cv_vals)]
  
  # GCV
  gcv_vals <- GCV_loclin(X, Y, h = h_grid, kernel = ker)
  h_gcv <- h_grid[which.min(gcv_vals)]
  
  res_list[[ker]] <- data.frame(
    Kernel = ker,
    h_RT   = h_rt,
    h_CV   = h_cv,
    h_GCV  = h_gcv
  )
}

result <- do.call(rbind, res_list)


#print(result)
#  Kernel      h_RT      h_CV     h_GCV
#1  gauss 0.1700174 0.1843434 0.1747475
#2   epan 0.2727964 0.4338384 0.2419192
h_RT_gauss  <- result$h_RT[result$Kernel == "gauss"]
h_CV_gauss  <- result$h_CV[result$Kernel == "gauss"]
h_GCV_gauss <- result$h_GCV[result$Kernel == "gauss"]

h_RT_epan  <- result$h_RT[result$Kernel == "epan"]
h_CV_epan  <- result$h_CV[result$Kernel == "epan"]
h_GCV_epan <- result$h_GCV[result$Kernel == "epan"]

#c
x_grid <- seq(min(X), max(X), length.out = 200)
m_gauss_RT  <- loclin_reg(X, Y, x_grid, h_RT_gauss, kernel = "gauss")
m_gauss_CV  <- loclin_reg(X, Y, x_grid, h_CV_gauss, kernel = "gauss")
m_gauss_GCV <- loclin_reg(X, Y, x_grid, h_GCV_gauss, kernel = "gauss")

m_epan_RT  <- loclin_reg(X, Y, x_grid, h_RT_epan, kernel = "epan")
m_epan_CV  <- loclin_reg(X, Y, x_grid, h_CV_epan, kernel = "epan")
m_epan_GCV <- loclin_reg(X, Y, x_grid, h_GCV_epan, kernel = "epan")

#plot(X, Y, pch = 16, col = "grey",
#main = "Local Linear Regression – Gaussian kernel",
#xlab = "X", ylab = "Y")

#lines(x_grid, m_gauss_RT,  col = "blue",  lwd = 2)
#lines(x_grid, m_gauss_CV,  col = "red",   lwd = 2)
#lines(x_grid, m_gauss_GCV, col = "green", lwd = 2)

#legend("topright",
#legend = c("RT", "CV", "GCV"),
#col = c("blue", "red", "green"),
#lwd = 2)


#plot(X, Y, pch = 16, col = "grey",
#main = "Local Linear Regression – Epan kernel",
#xlab = "X", ylab = "Y")

#lines(x_grid, m_epan_RT,  col = "blue",  lwd = 2)
#lines(x_grid, m_epan_CV,  col = "red",   lwd = 2)
#lines(x_grid, m_epan_GCV, col = "green", lwd = 2)

#legend("topright",
#legend = c("RT", "CV", "GCV"),
#col = c("blue", "red", "green"),
#lwd = 2)

#d)

res_d <- data.frame(
  iter = integer(600),
  kernel = character(600),
  method = character(600),
  h = numeric(600),
  stringsAsFactors = FALSE
)

idx <- 1

for (i in 1:100) {
  X <- runif(100, 0, pi)
  sd_eps <- sqrt(m(X)^2 / 64)
  eps <- rnorm(100, 0, sd_eps)
  Y <- m(X) + eps
  h_grid <- seq(0.05, 1, length.out = 100)
  
  for (ker in c("gauss", "epan")) {
    h_rt <- h_ll_RT(X, Y, ker)
    h_cv <- h_grid[which.min(CV_loclin(X, Y, h_grid, ker))]
    h_gcv <- h_grid[which.min(GCV_loclin(X, Y, h_grid, ker))]
    
    res_d[idx, ] <- list(i, ker, "RT", h_rt); idx <- idx + 1
    res_d[idx, ] <- list(i, ker, "CV", h_cv); idx <- idx + 1
    res_d[idx, ] <- list(i, ker, "GCV", h_gcv); idx <- idx + 1
  }
}
# Tạo thư mục lưu plot nếu chưa có
dir.create("./output/plots", recursive = TRUE, showWarnings = FALSE)

# Lưu boxplot Gaussian kernel
png(filename = "./output/plots/ex_5d_Gaussian.png", width = 800, height = 600)
boxplot(h ~ method,
        data = subset(res_d, kernel == "gauss"),
        main = "Bandwidth variability (Gaussian kernel)",
        ylab = "Optimal bandwidth h",
        col = c("lightblue", "pink", "lightgreen"))
dev.off()

# Lưu boxplot Epan kernel
png(filename = "./output/plots/ex_5d_Epan.png", width = 800, height = 600)
boxplot(h ~ method,
        data = subset(res_d, kernel == "epan"),
        main = "Bandwidth variability (Epan kernel)",
        ylab = "Optimal bandwidth h",
        col = c("lightblue", "pink", "lightgreen"))
dev.off()