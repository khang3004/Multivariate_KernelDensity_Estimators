source("R/multivariate_kde/data_transformation.R")
# Giả sử bạn đã có hàm mkde_basic(x_train, x_query, h) ở bài trước

#' MKDE with Whitening (The "Bread-to-Sphere" Method)
#'
#' @param x_train Matrix (n x d). Training data.
#' @param x_query Matrix (m x d). Query points.
#' @param bandwidth Numeric or Vector.
#'        If NULL, it auto-selects bandwidth ON THE WHITENED DATA using Scott's Rule.
#'
#' @return Vector of density estimates.
#' @export
mkde_whitened <- function(x_train, x_query, bandwidth = NULL) {
  d <- ncol(x_train)

  # --- 1. LEARN TRANSFORMATION ---
  params <- get_whitening_params(x_train)

  # --- 2. TRANSFORM DATA (X -> Z) ---
  z_train <- whiten_data(x_train, params)
  z_query <- whiten_data(x_query, params)

  # --- 3. BANDWIDTH SELECTION (On Z) ---
  # Vì Z đã tròn quay (variance ~ 1), ta có thể dùng Scott's Rule đơn giản
  if (is.null(bandwidth)) {
    n <- nrow(z_train)
    # Scott's Rule for Product Kernel: h = n^(-1/(d+4))
    # Vì Z đã chuẩn hóa variance=1, ta ko cần nhân với sigma nữa
    h_opt <- n^(-1/(d + 4))
    bandwidth <- rep(h_opt, d)
  }

  # --- 4. RUN CORE MKDE ON Z ---
  # Gọi hàm MKDE cơ bản mà bạn đã viết (Product Kernel)
  # Lưu ý: Hàm này tính density f_Z(z)
  # Bạn cần đảm bảo đã load hàm mkde_core cơ bản
  dens_z <- mkde_product_kernel(z_train, z_query, h = bandwidth)

  # --- 5. JACOBIAN CORRECTION (Back to X) ---
  # f_X(x) = f_Z(z) * det(W)
  # Thực ra công thức chuẩn là: f_X = f_Z * det(Jacobian of z->x)^-1
  # Z = (X-mu)W  => dZ/dX = W. => dX/dZ = W^-1.
  # Density transform: f_X(x) = f_Z(z(x)) * |det(W)|

  det_W <- det(params$W)
  dens_x <- dens_z * abs(det_W)

  return(dens_x)
}

# --- HÀM MKDE CƠ BẢN (Nhắc lại để code chạy được) ---
mkde_product_kernel <- function(x_train, x_query, h) {
  n <- nrow(x_train)
  m <- nrow(x_query)
  d <- ncol(x_train)

  # Nếu h là số, nhân bản lên d chiều
  if (length(h) == 1) h <- rep(h, d)

  # Product Kernel: Gaussian
  # K(u) = product( (1/sqrt(2pi)) * exp(-0.5 * u_j^2) )

  # Để tối ưu, ta tính log-density rồi exp lại để tránh underflow
  const <- -0.5 * log(2 * pi)

  preds <- numeric(m)

  for (i in 1:m) {
    # Tính khoảng cách đã scale: (x_train - query_i) / h
    # Matrix (n x d)
    diffs <- sweep(x_train, 2, x_query[i, ], "-")
    scaled <- sweep(diffs, 2, h, "/")

    # Tính log kernel cho từng chiều: -0.5 * u^2 - log(h) + const
    log_k <- -0.5 * scaled^2 + const
    # Trừ log(h) cho từng chiều (chuẩn hóa 1/h)
    log_k <- sweep(log_k, 2, log(h), "-")

    # Cộng tổng các chiều (Product Kernel -> Sum of Logs)
    log_K_multi <- rowSums(log_k)

    # Tính trung bình mật độ của n điểm training
    # Mean of exp(log_K)
    # Dùng log-sum-exp trick nếu cần, ở đây viết đơn giản:
    preds[i] <- mean(exp(log_K_multi))
  }

  return(preds)
}