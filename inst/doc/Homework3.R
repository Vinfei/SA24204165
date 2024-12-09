## -----------------------------------------------------------------------------
# 设置参数
set.seed(123)
n <- 1000       # 样本大小
num_sim <- 10000  # 蒙特卡洛模拟次数

# 定义函数计算样本偏度
sample_skewness <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sd(x)
  m3 <- sum((x - mean_x)^3) / n
  b1 <- (n / ((n - 1) * (n - 2))) * (m3 / (sd_x^3))
  return(b1)
}

# 初始化向量存储sqrt(b1)
sqrt_b1 <- numeric(num_sim)

# 蒙特卡洛模拟
for (i in 1:num_sim) {
  sample <- rnorm(n) # 生成标准正态分布样本
  b1 <- sample_skewness(sample)
  
  # 为避免负值导致复数，取绝对值
  if (b1 >= 0) {
    sqrt_b1[i] <- sqrt(b1)
  } else {
    sqrt_b1[i] <- NA  # 标记为缺失值
  }
}

# 移除缺失值
sqrt_b1 <- sqrt_b1[!is.na(sqrt_b1)]

# 计算所需分位数
quantiles_sim <- quantile(sqrt_b1, probs = c(0.025, 0.05, 0.95, 0.975))

# 大样本近似的分位数
# 由于 sqrt(b1) ~ N(0, sqrt(6/n))，计算理论分位数
# 注意：由于我们取了绝对值，理论分布应考虑对称性
# 因此，计算双侧分位数
sd_approx <- sqrt(6 / n)
quantiles_theory <- qnorm(c(0.025, 0.05, 0.95, 0.975), mean = 0, sd = sd_approx)

# 计算标准误差 using Var(x̂_q) = q(1 - q) / (n f(x_q)^2 )
# 这里的f(x_q)是正态分布的密度函数
# 注意：因为我们取了绝对值，密度应取双侧的

# 计算f(x_q)
f_xq <- dnorm(quantiles_sim, mean = 0, sd = sd_approx)

# 计算Var(x_q)和标准误差
var_xq <- (c(0.025, 0.05, 0.95, 0.975) * (1 - c(0.025, 0.05, 0.95, 0.975))) / (n * (f_xq)^2)
se_xq <- sqrt(var_xq)

# 比较模拟分位数与理论分位数
comparison <- data.frame(
  Quantile = c(0.025, 0.05, 0.95, 0.975),
  Simulated = quantiles_sim,
  Theoretical = quantiles_theory,
  SE = se_xq
)

print("分位数比较:")
print(comparison)


## -----------------------------------------------------------------------------
# 设置参数
set.seed(123)
n <- 1000  # 样本数量
n_sim <- 1000  # 模拟次数

# Step 1: 双变量正态分布
library(MASS)
mu <- c(0, 0)
Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # 相关系数为 0.5 的协方差矩阵

# 初始化 p 值向量
p_values_pearson_normal <- numeric(n_sim)
p_values_spearman_normal <- numeric(n_sim)
p_values_kendall_normal <- numeric(n_sim)

# 模拟双变量正态分布
for (i in 1:n_sim) {
  data <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X <- data[, 1]
  Y <- data[, 2]
  
  # 相关性检验
  p_values_pearson_normal[i] <- cor.test(X, Y, method = "pearson")$p.value
  p_values_spearman_normal[i] <- cor.test(X, Y, method = "spearman")$p.value
  p_values_kendall_normal[i] <- cor.test(X, Y, method = "kendall")$p.value
}

# 计算正态分布下的拒绝率 (alpha = 0.05)
power_pearson_normal <- mean(p_values_pearson_normal < 0.05)
power_spearman_normal <- mean(p_values_spearman_normal < 0.05)
power_kendall_normal <- mean(p_values_kendall_normal < 0.05)

cat("双变量正态分布时的检验力:\n")
cat("Pearson检验力: ", power_pearson_normal, "\n")
cat("Spearman检验力: ", power_spearman_normal, "\n")
cat("Kendall检验力: ", power_kendall_normal, "\n")

# Step 2: 替代分布 - 非正态分布 (例如, X 和 Y 存在非线性关系)
p_values_pearson_alt <- numeric(n_sim)
p_values_spearman_alt <- numeric(n_sim)
p_values_kendall_alt <- numeric(n_sim)

for (i in 1:n_sim) {
  X <- rnorm(n)
  Y <- X^2 + rnorm(n, sd = 0.1)  # X 和 Y 存在非线性关系
  
  # 相关性检验
  p_values_pearson_alt[i] <- cor.test(X, Y, method = "pearson")$p.value
  p_values_spearman_alt[i] <- cor.test(X, Y, method = "spearman")$p.value
  p_values_kendall_alt[i] <- cor.test(X, Y, method = "kendall")$p.value
}

# 计算非正态分布下的拒绝率 (alpha = 0.05)
power_pearson_alt <- mean(p_values_pearson_alt < 0.05)
power_spearman_alt <- mean(p_values_spearman_alt < 0.05)
power_kendall_alt <- mean(p_values_kendall_alt < 0.05)

cat("\n替代分布（非正态分布）时的检验力:\n")
cat("Pearson检验力: ", power_pearson_alt, "\n")
cat("Spearman检验力: ", power_spearman_alt, "\n")
cat("Kendall检验力: ", power_kendall_alt, "\n")

## -----------------------------------------------------------------------------
### 1. 假设检验问题
# H0: 两种方法效能无显著差异 (mu1 = mu2)
# H1: 两种方法效能有显著差异 (mu1 != mu2)
# 这里 mu1 和 mu2 分别是两种方法的效能均值。

### 2. 使用哪种检验
# 1. Z检验: 适用于样本量很大且已知总体方差的情况下。这里未提供总体方差信息，因此Z检验不适合。
# 2. 两样本t检验: 适用于比较两组独立样本的均值是否有显著差异，假设方差不相等或相等的情况下都可以用。
# 3. 配对t检验: 适用于相同对象在两种不同处理下的效能比较（即配对数据）。如果两个效能值是针对同一组实验对象产生的，则应该使用配对t检验。
# 4. McNemar检验: 适用于分类变量，尤其是配对数据的二项检验。这里效能是连续变量，所以不适用。
# 假设两种方法的效能是基于相同的实验对象测量的（即配对数据），配对t检验会是更适合的选择。
# 如果这两组效能是独立的（不同实验对象），则使用两样本t检验。

### 3. 进行假设检验所需的最少必要信息
# - 每种方法的效能值（样本均值、方差）
# - 样本大小（10,000次实验）
# - 显著性水平 (\(\alpha = 0.05\))
# - 是否配对数据


# 设置随机种子以便结果可重现
set.seed(123)

# 样本大小，即实验的次数
n <- 10000  

# 模拟两组效能数据:
# 方法1的效能均值为0.651，标准差为0.05
# 方法2的效能均值为0.676，标准差为0.05
method1 <- rnorm(n, mean = 0.651, sd = 0.05)  # 第一组方法的效能数据
method2 <- rnorm(n, mean = 0.676, sd = 0.05)  # 第二组方法的效能数据

# 1. 配对t检验（假设实验是基于相同对象）
# H0: 两种方法效能无显著差异 (mu1 = mu2)
# H1: 两种方法效能有显著差异 (mu1 != mu2)
paired_t_test_result <- t.test(method1, method2, paired = TRUE)

# 打印配对t检验结果
print("配对t检验结果：")
print(paired_t_test_result)

# 2. 两样本t检验（假设两组数据来自不同实验对象）
# H0: 两种方法效能无显著差异 (mu1 = mu2)
# H1: 两种方法效能有显著差异 (mu1 != mu2)
independent_t_test_result <- t.test(method1, method2, paired = FALSE)

# 打印两样本t检验结果
print("两样本t检验结果：")
print(independent_t_test_result)

# 结论：
# 比较两个t检验的p值，如果p值小于0.05，我们可以拒绝原假设，说明效能存在显著差异。

