set.seed(1)
# 5
# a
m <- function(x) {
  1 + sin((x^2)) / x^2
}

# X pp deu tren (0,pi) : Uniform(0,pi)
X <- runif(100, 0, pi)

# Độ lệch chuẩn của ei: sqrt(m(x)^2 / 64)
sd_eps <- sqrt(m(X)^2 / 64)

eps <- rnorm(100, 0, sd_eps)

Y <- m(X) + eps

# b
kernel_fun <- function(u, kernel) {
  if (kernel == "gauss") {
    return(exp(-u^2 / 2) / sqrt(2 * pi))
  } else if (kernel == "epan") {
    return(0.75 * (1 - u^2) * (abs(u) <= 1))
  } else {
    stop("Unknown kernel")
  }
}

loclin_reg <- function(x, y, x_eval, h, kernel) {
  u <- x - x_eval
  uh <- u / h
  w <- kernel_fun(uh, kernel)
  s0 <- mean(w)
  s1 <- mean(w * (x - x_eval))
  s2 <- mean(w * (x - x_eval)^2)
  ww <- ((s2 - s1 * (x - x_eval)) * w) / (s2 * s0 - s1^2)
  y_hat <- mean(ww * y)
  return(y_hat)
}

loclin_reg <- Vectorize(FUN = loclin_reg, vectorize.args = "x_eval")

h_ll_RT <- function(x, y, kernel) {
  if (kernel == "gauss") {
    R_k <- 0.5 / sqrt(pi)
    sigma2_k <- 1
  } else if (kernel == "epan") {
    R_k <- 3 / 5
    sigma2_k <- 1 / 5
  } else {
    stop("Unknown kernel")
  }
  # Quartic fit
  mod_Q <- lm(y ~ poly(x, raw = TRUE, degree = 4))
  # Estimates of unknown quantities
  int_sigma2_hat <- diff(range(x)) * sum(mod_Q$residuals^2) / mod_Q$df.residual
  m_2prime <- (2 * mod_Q$coefficients[3] + 6 * mod_Q$coefficients[4] * x +
    12 * mod_Q$coefficients[5] * x^2)^2
  a_hat <- mean(m_2prime)
  # h_RT
  h_RT <- ((R_k * int_sigma2_hat) / (a_hat * sigma2_k * length(x)))^(1 / 5)
  return(h_RT)
}

## Cross-validation
CV_loclin <- function(x, y, h, kernel) {
  n <- length(y)
  S_diag <- numeric(n)
  for (i in 1:n) {
    u <- x - x[i]
    w <- kernel_fun(u / h, kernel)
    a0 <- mean(w)
    a1 <- mean(w * u)
    a2 <- mean(w * u^2)
    S_diag[i] <- (a2 * w[i]) / (a2 * a0 - a1^2) / n
  }
  y_hat <- loclin_reg(x = x, y = y, x_eval = x, h = h, kernel = kernel)
  loocv <- mean(((y - y_hat) / (1 - S_diag))^2)
  return(loocv)
}

CV_loclin <- Vectorize(FUN = CV_loclin, vectorize.args = "h")


GCV_loclin <- function(x, y, h, kernel) {
  n <- length(y)
  S_diag <- numeric(n)
  for (i in 1:n) {
    u <- x - x[i]
    w <- kernel_fun(u / h, kernel)
    a0 <- mean(w)
    a1 <- mean(w * u)
    a2 <- mean(w * u^2)
    S_diag[i] <- (a2 * w[i]) / (a2 * a0 - a1^2) / n
  }
  y_hat <- loclin_reg(x = x, y = y, x_eval = x, h = h, kernel = kernel)
  gcv <- sum((y - y_hat)^2) / (n - sum(S_diag))^2
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

# print(result)
#  Kernel      h_RT      h_CV     h_GCV
# 1  gauss 0.1700174 0.1843434 0.1747475
# 2   epan 0.2727964 0.4338384 0.2419192
h_RT_gauss <- result$h_RT[result$Kernel == "gauss"]
h_CV_gauss <- result$h_CV[result$Kernel == "gauss"]
h_GCV_gauss <- result$h_GCV[result$Kernel == "gauss"]

h_RT_epan <- result$h_RT[result$Kernel == "epan"]
h_CV_epan <- result$h_CV[result$Kernel == "epan"]
h_GCV_epan <- result$h_GCV[result$Kernel == "epan"]

# c
x_grid <- seq(min(X), max(X), length.out = 200)
m_gauss_RT <- loclin_reg(X, Y, x_grid, h_RT_gauss, kernel = "gauss")
m_gauss_CV <- loclin_reg(X, Y, x_grid, h_CV_gauss, kernel = "gauss")
m_gauss_GCV <- loclin_reg(X, Y, x_grid, h_GCV_gauss, kernel = "gauss")

m_epan_RT <- loclin_reg(X, Y, x_grid, h_RT_epan, kernel = "epan")
m_epan_CV <- loclin_reg(X, Y, x_grid, h_CV_epan, kernel = "epan")
m_epan_GCV <- loclin_reg(X, Y, x_grid, h_GCV_epan, kernel = "epan")

# plot(X, Y, pch = 16, col = "grey",
# main = "Local Linear Regression – Gaussian kernel",
# xlab = "X", ylab = "Y")

# lines(x_grid, m_gauss_RT,  col = "blue",  lwd = 2)
# lines(x_grid, m_gauss_CV,  col = "red",   lwd = 2)
# lines(x_grid, m_gauss_GCV, col = "green", lwd = 2)

# legend("topright",
# legend = c("RT", "CV", "GCV"),
# col = c("blue", "red", "green"),
# lwd = 2)


# plot(X, Y, pch = 16, col = "grey",
# main = "Local Linear Regression – Epan kernel",
# xlab = "X", ylab = "Y")

# lines(x_grid, m_epan_RT,  col = "blue",  lwd = 2)
# lines(x_grid, m_epan_CV,  col = "red",   lwd = 2)
# lines(x_grid, m_epan_GCV, col = "green", lwd = 2)

# legend("topright",
# legend = c("RT", "CV", "GCV"),
# col = c("blue", "red", "green"),
# lwd = 2)


# d)

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
      data.frame(iter = b, kernel = ker, method = "RT", h = h_rt),
      data.frame(iter = b, kernel = ker, method = "CV", h = h_cv),
      data.frame(iter = b, kernel = ker, method = "GCV", h = h_gcv)
    )
  }
}
boxplot(h ~ method,
  data = subset(res_d, kernel == "gauss"),
  main = "Bandwidth variability (Gaussian kernel)",
  ylab = "Optimal bandwidth h",
  col = c("lightblue", "pink", "lightgreen")
)
boxplot(h ~ method,
  data = subset(res_d, kernel == "epan"),
  main = "Bandwidth variability (Epan kernel)",
  ylab = "Optimal bandwidth h",
  col = c("lightblue", "pink", "lightgreen")
)
