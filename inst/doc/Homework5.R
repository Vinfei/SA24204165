## -----------------------------------------------------------------------------
# 加载库
library(bootstrap)
attach(scor)

# 将数据转换为矩阵
x <- as.matrix(scor)

# 获取数据的行数
n <- nrow(x)

# 初始化Jackknife估计的向量
theta_jack <- numeric(n)

# 计算原始数据协方差矩阵的特征值
lambda <- eigen(cov(x))$values
theta_hat <- max(lambda) / sum(lambda)

# 对每个观测值进行Jackknife
for (i in 1:n) {
  # 排除第i个观测值
  y <- x[-i, ]
  # 计算协方差矩阵
  s <- cov(y)
  # 计算特征值
  lambda <- eigen(s)$values
  # 计算Jackknife估计
  theta_jack[i] <- max(lambda) / sum(lambda)
}

# 计算Jackknife估计的偏差和标准误差
bias_jack <- (n - 1) * (mean(theta_jack) - theta_hat)
se_jack <- sqrt((n - 1) / n * sum((theta_jack - mean(theta_jack))^2))

# 输出结果
cat("Jackknife估计的偏差:", bias_jack, "\n")
cat("Jackknife估计的标准误差:", se_jack, "\n")

# 将结果存储在列表中
result <- list(est = theta_hat, bias = bias_jack, se = se_jack)

# 分离数据环境
detach(scor)
detach(package:bootstrap)

# 打印结果列表
print(result)

## -----------------------------------------------------------------------------
# 加载必要的库
library(DAAG)
attach(ironslag)

# 定义化学成分的范围
a <- seq(10, 40, 0.1)

# 调整绘图参数和尺寸
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(1, 1, 2, 1))

# 拟合四个模型并绘制数据与拟合曲线
# 线性模型
model1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main = "Linear Model", pch = 16, xlab = "Chemical", ylab = "Magnetic")
lines(a, predict(model1, newdata = data.frame(chemical = a)), col = "blue")

# 二次多项式模型
model2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main = "Quadratic Model", pch = 16, xlab = "Chemical", ylab = "Magnetic")
lines(a, predict(model2, newdata = data.frame(chemical = a)), col = "red")

# 对数模型
model3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main = "Exponential Model", pch = 16, xlab = "Chemical", ylab = "Magnetic")
lines(a, exp(predict(model3, newdata = data.frame(chemical = a))), col = "green")

# 三次多项式模型
model4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))
plot(chemical, magnetic, main = "Cubic Model", pch = 16, xlab = "Chemical", ylab = "Magnetic")
lines(a, predict(model4, newdata = data.frame(chemical = a)), col = "purple")

# 恢复单一绘图模式
par(mfrow = c(1, 1))

# 计算每个模型的调整后的R²
rsq_values <- sapply(list(model1, model2, model3, model4), function(model) {
  summary(model)$adj.r.squared
})

# 输出调整后的R²值
print(rsq_values)

# 留一法交叉验证评估预测误差
n <- length(magnetic)
errors <- matrix(0, ncol = 4, nrow = n)

for (k in 1:n) {
  temp_data <- ironslag[-k, ]
  full_data <- ironslag[k, ]

  model1_fit <- lm(magnetic ~ chemical, data = temp_data)
  errors[k, 1] <- full_data$magnetic - predict(model1_fit, newdata = full_data)

  model2_fit <- lm(magnetic ~ chemical + I(chemical^2), data = temp_data)
  errors[k, 2] <- full_data$magnetic - predict(model2_fit, newdata = full_data)

  model3_fit <- lm(log(magnetic) ~ chemical, data = temp_data)
  errors[k, 3] <- full_data$magnetic - exp(predict(model3_fit, newdata = full_data))

  model4_fit <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3), data = temp_data)
  errors[k, 4] <- full_data$magnetic - predict(model4_fit, newdata = full_data)
}

# 计算每个模型的预测误差
prediction_errors <- colMeans(errors^2)
print(prediction_errors)

# 模型(2)具有最小的预测误差，因此二次模型(2)也通过交叉验证被选中。

# 分离数据环境
detach(ironslag)
detach(package:DAAG)


## -----------------------------------------------------------------------------
# 定义函数
cvm_test <- function(sample1, sample2, num_permutations = 199) {
  len1 <- length(sample1)
  len2 <- length(sample2)
  combined <- c(sample1, sample2)
  total_len <- len1 + len2
  
  # 计算经验分布函数
  ecdf1 <- numeric(total_len)
  ecdf2 <- numeric(total_len)
  
  for (i in 1:total_len) {
    ecdf1[i] <- mean(as.integer(combined[i] <= sample1))
    ecdf2[i] <- mean(as.integer(combined[i] <= sample2))
  }
  
  # 计算观察到的统计量
  observed_statistic <- ((len1 * len2) / total_len^2) * sum((ecdf1 - ecdf2)^2)
  
  # 进行置换检验
  permuted_statistics <- replicate(num_permutations, {
    shuffled_indices <- sample(1:total_len)
    shuffled_combined <- combined[shuffled_indices]
    sample1_permuted <- shuffled_combined[1:len1]
    sample2_permuted <- shuffled_combined[(len1 + 1):total_len]
    
    # 重新计算经验分布函数
    for (i in 1:total_len) {
      ecdf1[i] <- mean(as.integer(shuffled_combined[i] <= sample1_permuted))
      ecdf2[i] <- mean(as.integer(shuffled_combined[i] <= sample2_permuted))
    }
    
    # 计算置换后的统计量
    ((len1 * len2) / total_len^2) * sum((ecdf1 - ecdf2)^2)
  })
  
  # 计算p值
  permuted_statistics <- c(permuted_statistics, observed_statistic)
  p_value <- mean(permuted_statistics >= observed_statistic)
  
  return(list(statistic = observed_statistic, p_value = p_value))
}

# 加载数据
attach(chickwts)

# 提取数据
x1 <- weight[feed == "soybean"]
x2 <- weight[feed == "sunflower"]
x3 <- weight[feed == "linseed"]

# 运行检验
test_result1 <- cvm_test(x1, x3)
test_result2 <- cvm_test(x2, x3)

# 输出结果
print(test_result1)
print(test_result2)

# 比较CvM检验的p值不显著。无法表明这些分布之间存在差异。

# 卸载数据
detach(chickwts)

## -----------------------------------------------------------------------------
# 定义函数
spearman_permutation_test <- function(data_x, data_y) {
  # 计算相关系数
  original_test <- cor.test(data_x, data_y, method = "spearman")
  
  # 获取样本长度
  sample_size <- length(data_x)
  
  # 进行num_permutations次置换，计算每次置换后的斯皮尔曼秩相关系数
  permuted_results <- replicate(num_permutations, expr = {
    shuffled_index <- sample(1:sample_size)
    cor.test(data_x, data_y[shuffled_index], method = "spearman")$estimate
  })
  
  # 将原始统计量与置换统计量合并
  all_results <- c(original_test$estimate, permuted_results)
  
  # 计算p值，即原始统计量小于等于置换统计量的比例
  p_value <- mean(all_results >= original_test$estimate)
  
  # 返回结果
  return(list(rho = original_test$estimate, p_value = p_value))
}

# 加载包
library(MASS)

# 设置均值和协方差矩阵
mean_vector <- c(0, 0)
covariance_matrix <- matrix(c(1, 0.5, 0.5, 1), 2, 2)

# 设置样本大小和置换次数
sample_count <- 30
num_permutations <- 1000

# 生成双变量正态分布样本
normal_samples <- mvrnorm(sample_count, mean_vector, covariance_matrix)

# 应用斯皮尔曼秩相关检验和置换检验
normal_cor_test <- cor.test(normal_samples[, 1], normal_samples[, 2], method = "spearman")
normal_spearman_test <- spearman_permutation_test(normal_samples[, 1], normal_samples[, 2])

# 输出结果
print(normal_cor_test)
print(normal_spearman_test)

# 生成对数正态分布样本
lognormal_samples <- exp(mvrnorm(sample_count, mean_vector, covariance_matrix))

# 再次应用斯皮尔曼秩相关检验和置换检验
lognormal_cor_test <- cor.test(lognormal_samples[, 1], lognormal_samples[, 2], method = "spearman")
lognormal_spearman_test <- spearman_permutation_test(lognormal_samples[, 1], lognormal_samples[, 2])

# 输出结果
print(lognormal_cor_test)
# 两个测试的p值都是显著的，并且非常接近。
print(lognormal_spearman_test)
# 两个测试的p值都是显著的，并且非常接近。

# 卸载MASS包
detach(package:MASS)

