## -----------------------------------------------------------------------------
# 加载库
library(ggplot2)

# 定义生成瑞利分布样本的函数
generate_rayleigh_samples <- function(n, sigma) {
  U <- runif(n)  
  X <- sigma * sqrt(-2 * log(1 - U))  # 逆变换抽样法生成瑞利分布样本
  return(X)
}

# 定义函数计算众数（使用直方图的最高峰作为众数的近似值）
estimate_mode <- function(samples) {
  dens <- density(samples)  # 样本密度估计
  mode_value <- dens$x[which.max(dens$y)]  # 找到密度最高点
  return(mode_value)
}

# 生成样本并验证众数是否接近理论众数
check_mode_for_sigma <- function(sigma, n = 10000) {
  # 生成瑞利分布样本
  samples <- generate_rayleigh_samples(n, sigma)
  
  # 计算理论众数
  theoretical_mode <- sigma
  
  # 估计样本众数
  estimated_mode <- estimate_mode(samples)
  
  # 绘制直方图与理论众数对比
  hist(samples, breaks = 50, probability = TRUE, 
       main = paste("Rayleigh Distribution (σ =", sigma, ")"),
       xlab = "Value", col = "lightblue", border = "white")
  
  # 添加样本密度曲线
  lines(density(samples), col = "blue", lwd = 2)
  
  # 添加理论众数和估计众数的标记
  abline(v = theoretical_mode, col = "red", lwd = 2, lty = 2)
  abline(v = estimated_mode, col = "green", lwd = 2, lty = 2)
  
  # 显示理论众数和估计众数
  legend("topright", legend = c(paste("Theoretical Mode =", round(theoretical_mode, 2)),
                                paste("Estimated Mode =", round(estimated_mode, 2))),
         col = c("red", "green"), lty = 2, lwd = 2)
  
  # 返回估计众数与理论众数
  return(list(theoretical_mode = theoretical_mode, estimated_mode = estimated_mode))
}

# 设定不同的 σ 值进行验证
sigmas <- c(1, 2, 3)

# 对每个 σ 值进行检查
for (sigma in sigmas) {
  result <- check_mode_for_sigma(sigma)
  cat("For σ =", sigma, ": Theoretical Mode =", result$theoretical_mode, 
      "Estimated Mode =", result$estimated_mode, "\n")
}

## -----------------------------------------------------------------------------
library(ggplot2)

# 定义生成混合正态分布样本的函数
generate_mixture_samples <- function(n, p1, mean1 = 0, sd1 = 1, mean2 = 3, sd2 = 1) {
  # 生成随机分布的成分选择
  component <- rbinom(n, 1, p1)  # 生成 0 或 1，概率 p1 来自 N(0,1)，概率 p2 来自 N(3,1)
  
  # 根据成分选择生成相应的正态分布样本
  samples <- ifelse(component == 1, rnorm(n, mean1, sd1), rnorm(n, mean2, sd2))
  
  return(samples)
}

# 定义函数绘制样本直方图并叠加密度曲线
plot_mixture_distribution <- function(samples, p1) {
  df <- data.frame(samples = samples)
  
  # 绘制直方图和密度曲线
  ggplot(df, aes(x = samples)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
    geom_density(color = "red", size = 1) +
    labs(title = paste("Mixture of N(0,1) and N(3,1) with p1 =", p1),
         x = "Sample Value", y = "Density") +
    theme_minimal()
}

# 设定样本与不同的 p1 值
n <- 1000
p1_values <- c(0.2, 0.4,0.6, 0.75)

# 生成混合分布样本并绘制直方图
for (p1 in p1_values) {
  samples <- generate_mixture_samples(n, p1)
  print(plot_mixture_distribution(samples, p1))
}

## -----------------------------------------------------------------------------
# 设置参数
lambda <- 2       # 泊松过程的参数 λ
t <- 10           # 时间 t
shape <- 3        # Gamma 分布的 shape 参数
rate <- 1         # Gamma 分布的 rate 参数 (通常是 1/beta)

# 模拟泊松过程 N(t)
N_t <- rpois(10000, lambda * t)  # 生成10000个泊松随机数作为样本

# 模拟 Gamma 分布 Y_i
Y_i <- function(n) {
  rgamma(n, shape = shape, rate = rate)
}

# 计算复合泊松过程 X(t)
X_t <- sapply(N_t, function(n) sum(Y_i(n)))  # 通过泊松过程次数求和 Gamma 随机变量

# 计算 X(t) 的经验均值和方差
mean_X_t <- mean(X_t)
var_X_t <- var(X_t)

# 理论均值和方差
E_Y <- shape / rate  # Gamma 分布的期望
Var_Y <- shape / rate^2  # Gamma 分布的方差
theoretical_mean_X_t <- lambda * t * E_Y
theoretical_var_X_t <- lambda * t * (Var_Y + E_Y^2)

# 输出结果
cat("模拟结果：\n")
cat("X(10) 的均值 =", mean_X_t, "\n")
cat("X(10) 的方差 =", var_X_t, "\n\n")

cat("理论值：\n")
cat("X(10) 的理论均值 =", theoretical_mean_X_t, "\n")
cat("X(10) 的理论方差 =", theoretical_var_X_t, "\n")

