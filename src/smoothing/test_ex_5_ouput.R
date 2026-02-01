source('src/smoothing/Bai_5.R')

n <- 100
x <- seq(0, 10, length.out = n)
y <- rep(100, n) # Y luôn bằng 100

# Chọn x_eval tại điểm giữa, h tùy ý
pred <- loclin_reg(x = x, y = y, x_eval = 5, h = 1, kernel = "gauss")

# 4. Kiểm tra
print(paste("Giá trị thực tế:", 100))
print(paste("Giá trị mem dự đoán:", pred))

if (abs(pred - 100) < 1e-5) {
  print("✅ MEM ĐÚNG! (Phép toán đã tự triệt tiêu)")
} else {
  print(paste("❌ MEM SAI! Dự đoán lệch. Tỷ lệ sai số:", 100/pred))
}