source("local_linear_regression.R")
source("bandwidth_smoothing_selection/cross_validation.R")
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
h_ll_RT <- function(x, y, kernel = c("gauss", "epan")) {
  kernel <- match.arg(kernel)
  
  if (kernel == "gauss") {
    R_k <- 0.5 / sqrt(pi)
    sigma2_k <- 1
  } else if (kernel == "epan") {
    R_k <- 3 / 5
    sigma2_k <- 1 / 5
  }
  
  # Quartic polynomial fit
  mod_Q <- lm(y ~ poly(x, raw = TRUE, degree = 4))
  
  # Estimate integrated variance
  int_sigma2_hat <- diff(range(x)) * 
    sum(residuals(mod_Q)^2) / mod_Q$df.residual
  
  # Estimate squared second derivative
  m_2prime_sq <- (
    2 * mod_Q$coefficients[3] +
      6 * mod_Q$coefficients[4] * x +
      12 * mod_Q$coefficients[5] * x^2
  )^2
  
  a_hat <- mean(m_2prime_sq)
  
  # Rule-of-Thumb bandwidth
  h_RT <- ((R_k * int_sigma2_hat) /
             (a_hat * sigma2_k * length(x)))^(1/5)
  
  return(h_RT)
}

# Rule-of-Thumb bandwidths
h_rt_epan  <- h_ll_RT(X, Y, kernel = "epan")
h_rt_gauss <- h_ll_RT(X, Y, kernel = "gauss")
h_grid <- seq(0.05, 1.2, length.out = 50)

# Epanechnikov with function CV and GCV shortcut
cv_epan  <- cv_shortcut(X, Y, h_grid, kernel = "epan")
gcv_epan <- generalized_cv_shortcut(X, Y, h_grid, kernel = "epan")

# Gaussian with function CV and GCV shortcut
cv_gauss  <- cv_shortcut(X, Y, h_grid, kernel = "gauss")
gcv_gauss <- generalized_cv_shortcut(X, Y, h_grid, kernel = "gauss")

#create Data to restore and input
data <- data.frame(
  Kernel = rep(c("Epanechnikov", "Gaussian"), each = 3),
  Method = rep(c("Rule-of-Thumb", "CV", "GCV"), times = 2),
  h_opt = c(
    h_rt_epan,
    cv_epan$optimal_h,
    gcv_epan$optimal_h,
    h_rt_gauss,
    cv_gauss$optimal_h,
    gcv_gauss$optimal_h
  )
)
#   Kernel          Rule.of.Thumb  CV       GCV
#1 Epanechnikov     0.2727964 0.4489796 0.2377551
#2     Gaussian     0.1700174 0.1908163 0.1673469

#c)

x_grid <- seq(0, pi, length.out = 200)

# Local Linear fits for Epanechnikov
fit_epan_rt  <- local_linear_regression(X, Y, x_grid, h = h_rt_epan, kernel = "epan")
fit_epan_cv  <- local_linear_regression(X, Y, x_grid, h = cv_epan$optimal_h, kernel = "epan")
fit_epan_gcv <- local_linear_regression(X, Y, x_grid, h = gcv_epan$optimal_h, kernel = "epan")

# Local Linear fits for Gaussian
fit_gauss_rt  <- local_linear_regression(X, Y, x_grid, h = h_rt_gauss, kernel = "gauss")
fit_gauss_cv  <- local_linear_regression(X, Y, x_grid, h = cv_gauss$optimal_h, kernel = "gauss")
fit_gauss_gcv <- local_linear_regression(X, Y, x_grid, h = gcv_gauss$optimal_h, kernel = "gauss")


# Tạo thư mục lưu hình nếu chưa có
dir.create("./output/plots", recursive = TRUE, showWarnings = FALSE)

# Lưu 2 plot vào 1 file PNG, chia 1 hàng 2 cột
png(filename = "./output/plots/ex_5c_kernels.png", width = 1200, height = 600)

par(mfrow = c(1, 2))  # 1 row, 2 columns

# --- Plot 1: Gaussian Kernel ---
plot(X, Y, pch = 16, col = "grey",
     main = "Gaussian Kernel - Local Linear Regression",
     xlab = "X", ylab = "Y")
lines(x_grid, m(x_grid), col = "black", lwd = 2)           # True m(x)
lines(x_grid, fit_gauss_rt, col = "blue", lwd = 2)         # RT
lines(x_grid, fit_gauss_cv, col = "green", lwd = 2)        # CV
lines(x_grid, fit_gauss_gcv, col = "pink", lwd = 2)        # GCV
legend("topright",
       legend = c("True m(x)", "RT", "CV", "GCV"),
       col = c("black", "blue", "green", "pink"),
       lwd = 2,
       lty = 1)  # tất cả nét liền

# --- Plot 2: Epanechnikov Kernel ---
plot(X, Y, pch = 16, col = "grey",
     main = "Epanechnikov Kernel - Local Linear Regression",
     xlab = "X", ylab = "Y")
lines(x_grid, m(x_grid), col = "black", lwd = 2)           # True m(x)
lines(x_grid, fit_epan_rt, col = "blue", lwd = 2)          # RT
lines(x_grid, fit_epan_cv, col = "green", lwd = 2)         # CV
lines(x_grid, fit_epan_gcv, col = "pink", lwd = 2)         # GCV
legend("topright",
       legend = c("True m(x)", "RT", "CV", "GCV"),
       col = c("black", "blue", "green", "pink"),
       lwd = 2,
       lty = 1)  # tất cả nét liền

dev.off()  # đóng file PNG và lưu
#d)Repeat 100 times

# Create vectors to store bandwidths
h_epan_rt  <- numeric(100)
h_epan_cv  <- numeric(100)
h_epan_gcv <- numeric(100)

h_gauss_rt  <- numeric(100)
h_gauss_cv  <- numeric(100)
h_gauss_gcv <- numeric(100)

for (i in 1:100) {
  # Generate new dataset
  Xi <- runif(100, 0, pi)
  sd_eps <- sqrt(m(Xi)^2 / 64)
  Yi <- m(Xi) + rnorm(100, 0, sd_eps)
  
  # Rule-of-Thumb
  h_epan_rt[i]  <- h_ll_RT(Xi, Yi, "epan")
  h_gauss_rt[i] <- h_ll_RT(Xi, Yi, "gauss")
  
  # Cross-Validation
  h_epan_cv[i]  <- cv_shortcut(Xi, Yi, h_grid, "epan")$optimal_h
  h_gauss_cv[i] <- cv_shortcut(Xi, Yi, h_grid, "gauss")$optimal_h
  
  # Generalized Cross-Validation
  h_epan_gcv[i]  <- generalized_cv_shortcut(Xi, Yi, h_grid, "epan")$optimal_h
  h_gauss_gcv[i] <- generalized_cv_shortcut(Xi, Yi, h_grid, "gauss")$optimal_h
}

# Boxplot for Epanechnikov kernel
png(filename = "./output/plots/ex_5d_boxplots.png", width = 1200, height = 600)

# Chia canvas: 1 hàng 2 cột
par(mfrow = c(1, 2))

# --- Epanechnikov Kernel ---
boxplot(
  h_epan_rt, h_epan_cv, h_epan_gcv,
  names = c("RT", "CV", "GCV"),
  ylab = "Optimal bandwidth h",
  main = "Variability of Epanechnikov Bandwidths",
  col = c("blue", "green", "pink")
)

# --- Gaussian Kernel ---
boxplot(
  h_gauss_rt, h_gauss_cv, h_gauss_gcv,
  names = c("RT", "CV", "GCV"),
  ylab = "Optimal bandwidth h",
  main = "Variability of Gaussian Bandwidths",
  col = c("blue", "green", "pink")
)

dev.off()