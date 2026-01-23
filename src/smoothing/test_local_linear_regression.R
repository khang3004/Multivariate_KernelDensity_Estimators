# test_comparison.R
source("src/smoothing/kernel_regression.R")      # Cái cũ
source("src/smoothing/local_linear_regression.R") # Cái mới

# Load Data
data(mcycle, package = "MASS")
x <- mcycle$times
y <- mcycle$accel

# Grid
x_grid <- seq(min(x), max(x), length.out = 200)

# Run Models (Cùng bandwidth để công bằng)
h_val <- 2.0
y_nw <- nadaraya_watson(x, y, x_grid, h = h_val)
y_ll <- local_linear_regression(x, y, x_grid, h = h_val)

# Visualization
pdf("Compare_NW_vs_LocalLinear.pdf", width = 10, height = 6)
plot(x, y, pch=16, col=adjustcolor("gray", 0.5),
     main="Boundary Effect: NW vs Local Linear", xlab="Time", ylab="Accel")

# Vẽ Nadaraya-Watson (Màu xanh - Bị lệch ở biên)
lines(x_grid, y_nw, col="blue", lwd=2, lty=2)

# Vẽ Local Linear (Màu đỏ - Bám sát biên)
lines(x_grid, y_ll, col="red", lwd=2)

legend("topleft", legend=c("Nadaraya-Watson", "Local Linear"),
       col=c("blue", "red"), lty=c(2, 1), lwd=2)
dev.off()