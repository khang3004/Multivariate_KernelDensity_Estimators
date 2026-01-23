#load modules
source('R/multivariate_kde/bw_selection_multivariate.R')
source('R/multivariate_kde/mkde_core.R')



# 1. Simulate 2D Data (Bimodal Distribution)
set.seed(2026)
n <- 500

#Mode 1: Centered at (2,2)
data1 <- matrix(rnorm(n, mean = 2, sd = 1), ncol = 2)
head(data1)
#Mode 2: Centered at (-2,2)
data2 <- matrix(rnorm(n, mean = -2, sd = 1), ncol = 2)

#Combine
data_train <- rbind(data1, data2)
head(data_train)
colnames(data_train) <- c("X", "Y")

# 2. Bandwidth Selection (Scott's Rule or Silverman's Rule)
h_opt <- multi_silverman_bandwidth(data_train)
print("Optimal Bandwidths (Scott's Rule")
print(h_opt)

# Create Grid for Visualization
grid_size <- 50
x_seq <- seq(min(data_train[,1])-1, max(data_train[,1])+1,length.out=grid_size)
y_seq <- seq(min(data_train[,2])-1, max(data_train[,2])+1, length.out=grid_size)

# Create meshgrid (all combinations of x & y)
grid_matrix <- expand.grid(X=x_seq, Y=y_seq)
head(grid_matrix)
#4. Run MKDE
z_vals <- product_kernel_mkde(data = data_train, x_eval = as.matrix(grid_matrix), h_opt)
z_vals
# Reshape Z for contour plotting
z_matrix <- matrix(z_vals, nrow=grid_size, ncol=grid_size)
head(z_matrix)
#5. Visualization (Contour & Perspective)
pdf('report/figures/Multivariate_KDE_Report.pdf', width=10, height=5)
par(mfrow=c(1,2))

# Plot 1: Contour Plot with Data Points
image(x_seq, y_seq, z_matrix, col = terrain.colors(20),
      main = "2D KDE Contour (Product Kernel)", xlab = "X", ylab = "Y")
points(data_train, pch = 16, cex = 0.5, col = adjustcolor("black", 0.3))
contour(x_seq, y_seq, z_matrix, add = TRUE)

# Plot 2: 3D Surface
persp(x_seq, y_seq, z_matrix, theta = 30, phi = 30,
      col = "lightblue", shade = 0.5, border = NA,
      main = "3D Density Surface", zlab = "Density")

dev.off()
message("Report generated: Multivariate_KDE_Report.pdf")