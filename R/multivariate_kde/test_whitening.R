source("R/multivariate_kde/data_transformation.R")
source("R/multivariate_kde/mkde_whitened.R")
library(MASS) # Để sinh dữ liệu Multivariate Normal

# 1. Tạo dữ liệu "Bánh Mì" (Correlated Data)
set.seed(42)
Sigma <- matrix(c(1, 0.9, 0.9, 1), 2, 2) # Tương quan rất mạnh (0.9)
mu <- c(2, 5)
n <- 500
X <- mvrnorm(n, mu, Sigma)

# 2. Làm trắng (Whitening)
params <- get_whitening_params(X)
Z <- whiten_data(X, params)

# 3. Vẽ hình so sánh
par(mfrow=c(1, 2))
plot(X, main="Original (Bread Loaf)", col=rgb(0,0,1,0.5), pch=16, asp=1)
grid()
plot(Z, main="Whitened (Sphere)", col=rgb(1,0,0,0.5), pch=16, asp=1)
grid()





# --- TIẾP TỤC TỪ ĐOẠN CODE CỦA BẠN ---

# 1. Tạo lưới điểm (Grid) để vẽ contour
x_seq <- seq(min(X[,1])-1, max(X[,1])+1, length.out=50)
y_seq <- seq(min(X[,2])-1, max(X[,2])+1, length.out=50)
grid_pts <- expand.grid(x=x_seq, y=y_seq)
matrix_grid <- as.matrix(grid_pts)

# 2. CHẠY THỬ NGHIỆM: Naive vs Whitened

# Cách 1: Naive Product Kernel (Không Whitening)
# Chọn bừa bandwidth theo Scott's rule cho từng chiều độc lập
h_naive <- c(bandwidth.nrd(X[,1]), bandwidth.nrd(X[,2]))
dens_naive <- mkde_product_kernel(X, matrix_grid, h=h_naive)
z_naive <- matrix(dens_naive, nrow=50) # Reshape để vẽ

# Cách 2: Whitened MKDE (Cách xịn của bạn)
# Gọi hàm mkde_whitened bạn vừa viết (nó tự tính h bên trong)
dens_white <- mkde_whitened(X, matrix_grid)
z_white <- matrix(dens_white, nrow=50)

# 3. VẼ HÌNH SO SÁNH (THE REAL MAGIC)
par(mfrow=c(1, 2))

# Hình 1: Naive Kernel (Thất bại)
plot(X, main="Naive Product Kernel (Bad Fit)", col=rgb(0,0,1,0.2), pch=16, asp=1,
     xlab="X1", ylab="X2")
contour(x_seq, y_seq, z_naive, add=TRUE, col="red", lwd=2, nlevels=5)
# -> Nhận xét: Các vòng contour bị méo, nằm ngang dọc, không ôm sát dữ liệu nghiêng.

# Hình 2: Whitened MKDE (Hoàn hảo)
plot(X, main="Whitened MKDE (Perfect Fit)", col=rgb(0,0,1,0.2), pch=16, asp=1,
     xlab="X1", ylab="X2")
contour(x_seq, y_seq, z_white, add=TRUE, col="darkgreen", lwd=2, nlevels=5)
# -> Nhận xét: Các vòng contour nghiêng 45 độ, ôm khít lấy "ổ bánh mì".