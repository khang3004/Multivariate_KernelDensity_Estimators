# Load the refactored module
source('src/smoothing/kernel_regression.R')

# box::use(
#    src/smoothing/kernel_regresion
# )

#load mcycle data ---
data(mcycle, package = "MASS")
length(mcycle$accel)
head(mcycle$accel)

# 1. Setup Data & Grid
# Sort data for cleaner visulization (optional but good Practice)
ord <- order(mcycle$times)

x_data <- mcycle$times[ord]

y_data <- mcycle$accel[ord]
plot(x = x_data, y = y_data, pch = 16)
# Defin query grid
x_grid <- seq(min(x_data), max(x_data), length.out = 200)

# 2. Run Model
# Try bandwidth <- 2.3
y_pre <- nadaraya_watson(x_train = x_data,
                         y_train = y_data,
                         x_query = x_grid,
                         h = 2.3)


# 3. Visualization
# pdf("Nadaraya_Watson_Test.pdf", width = 10, height = 6)
# Base Plot (Raw Data)
graphics.off()
plot(x = x_data, y = y_data, # Add legend
     pch = 16, col = adjustcolor('gray',0.6),# Use transparency if supported, else just "gray"legend("topleft", legend = c('Raw Data', 'NW Estimator (h=2.3)'),
     main = "Nadaraya-Watson Regression on mcycl",
     xlab = "Time (ms)", ylab = "Head Acceleration (g)")

lines(x_grid, y_pre, col = 'red', lwd = 3)

# Add Legend
legend("topleft", legend = c("Raw Data", "NW Estimator (h=2.3)"),
       col = c("gray", "red"), pch = c(16, NA), lwd = c(NA, 3))

dev.off()
message('SUCCESS: Chart saved to src/smoothing/Nadaraya_Watson_Test.pdf')
































# In the lecture
u1 <- kernel_reg(x = mcycle$times, y = mcycle$accel, x_eval = 1, h = 0.5)

points(x = 30, y = u1, col = "blue")

x_plot <- seq(0, 60, length.out = 201)
y_hat <- kernel_reg(x = mcycle$times, y = mcycle$accel, x_eval = x_plot,
                    h = 0.5)

plot(x = mcycle, y = mcycle$accel, pch = 16, cex = 0.7)

lines(x = x_plot, y = y_hat, col = "blue")













