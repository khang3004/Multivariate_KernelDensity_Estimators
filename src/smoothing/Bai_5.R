set.seed(1)
#5
#a
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
kernel_fun <- function(u, kernel){
  if (kernel == "gauss") {
    return(exp(-u^2/2)/sqrt(2*pi))
  }
  else if (kernel == "epan") {
    return(0.75*(1 - u^2)*(abs(u) <= 1))
  }
  else stop("Unknown kernel")
}

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

## Cross-validation
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

result <- data.frame(
  Kernel = character(),
  h_RT   = numeric(),
  h_CV   = numeric(),
  h_GCV  = numeric()
)

for (ker in kernels) {
  
  # Rule of Thumb
  h_rt <- h_ll_RT(X, Y, kernel = ker)
  
  # CV
  cv_vals <- CV_loclin(X, Y, h = h_grid, kernel = ker)
  h_cv <- h_grid[which.min(cv_vals)]
  
  # GCV
  gcv_vals <- GCV_loclin(X, Y, h = h_grid, kernel = ker)
  h_gcv <- h_grid[which.min(gcv_vals)]
  
  result <- rbind(
    result,
    data.frame(
      Kernel = ker,
      h_RT   = h_rt,
      h_CV   = h_cv,
      h_GCV  = h_gcv
    )
  )
}

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
  iter   = integer(),
  kernel = character(),
  method = character(),
  h      = numeric()
)

for (i in 1:100) {
  
  # (a) Sinh dữ liệu mới
  X <- runif(100, 0, pi)
  sd_eps <- sqrt(m(X)^2 / 64)
  eps <- rnorm(100, 0, sd_eps)
  Y <- m(X) + eps
  
  h_grid <- seq(0.05, 1, length.out = 100)
  
  for (ker in c("gauss", "epan")) {
    
    # Rule of Thumb
    h_rt <- h_ll_RT(X, Y, ker)
    
    # CV
    cv_vals <- CV_loclin(X, Y, h_grid, ker)
    h_cv <- h_grid[which.min(cv_vals)]
    
    # GCV
    gcv_vals <- GCV_loclin(X, Y, h_grid, ker)
    h_gcv <- h_grid[which.min(gcv_vals)]
    
    # Lưu kết quả
    res_d <- rbind(
      res_d,
      data.frame(iter = i, kernel = ker, method = "RT",  h = h_rt),
      data.frame(iter = i, kernel = ker, method = "CV",  h = h_cv),
      data.frame(iter = i, kernel = ker, method = "GCV", h = h_gcv)
    )
  }
}
boxplot(h ~ method,
        data = subset(res_d, kernel == "gauss"),
        main = "Bandwidth variability (Gaussian kernel)",
        ylab = "Optimal bandwidth h",
        col = c("lightblue", "pink", "lightgreen"))
boxplot(h ~ method,
        data = subset(res_d, kernel == "epan"),
        main = "Bandwidth variability (Epan kernel)",
        ylab = "Optimal bandwidth h",
        col = c("lightblue", "pink", "lightgreen"))

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
