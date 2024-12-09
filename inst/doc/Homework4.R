## -----------------------------------------------------------------------------
# 设置参数
set.seed(123)  
N <- 1000  # 总假设数
N0 <- 950  # 原假设数
N1 <- 50   # 备择假设数
alpha <- 0.1
m <- 10000  

# 用于模拟生成 p 值的函数
generate_pvalues <- function(N0, N1) {
  p_values <- numeric(N)
  p_values[1:N0] <- runif(N0)
  p_values[(N0+1):N] <- rbeta(N1, 0.1, 1)
  return(p_values)
}

# Bonferroni 调整函数
bonferroni_adjust <- function(p_values, alpha) {
  p_adjusted <- p.adjust(p_values, method = "bonferroni")
  rejected <- p_adjusted <= alpha
  return(rejected)
}

# B-H 调整函数
bh_adjust <- function(p_values, alpha) {
  p_adjusted <- p.adjust(p_values, method = "BH")
  rejected <- p_adjusted <= alpha
  return(rejected)
}

# 进行模拟
simulations <- replicate(m, {
  p_values <- generate_pvalues(N0, N1)
  
  # 计算拒绝情况
  bonf_rejected <- bonferroni_adjust(p_values, alpha)
  bh_rejected <- bh_adjust(p_values, alpha)
  
  # 计算 FWER
  fwer_bonf <- mean(bonf_rejected[1:N0])  # 在原假设下的拒绝比例
  fwer_bh <- mean(bh_rejected[1:N0])
  
  # 计算 FDR
  fdr_bonf <- mean(bonf_rejected[(N0+1):N])  # 在备择假设下的拒绝比例
  fdr_bh <- mean(bh_rejected[(N0+1):N])
  
  # 计算 TPR
  tpr_bonf <- mean(bonf_rejected[(N0+1):N]) / N1  # 真正例率
  tpr_bh <- mean(bh_rejected[(N0+1):N]) / N1
  
  c(fwer_bonf, fwer_bh, fdr_bonf, fdr_bh, tpr_bonf, tpr_bh)
})

# 计算平均值
results <- apply(simulations, 1, mean)

# 输出结果
result_matrix <- matrix(results, nrow = 3, byrow = TRUE,
                         dimnames = list(c("FWER", "FDR", "TPR"),
                                         c("Bonferroni 校正", "B-H 校正")))
result_matrix


## -----------------------------------------------------------------------------
# 加载包
library(boot)

# 故障间隔时间数据
x <- aircondit[1]

# 定义最大似然估计函数
rate <- function(x, i) {
  # 从数据中提取 bootstrap 样本
  sample_data <- as.matrix(x[i, ])
  # 计算 MLE
  mle <- 1 / mean(sample_data)
  return(mle)
}

# 执行 bootstrap 过程
set.seed(123) # 设置随机种子以确保结果可重复
bootstrap_results <- boot(x, statistic = rate, R = 2000)

# 输出 bootstrap 统计结果
print(bootstrap_results)

## -----------------------------------------------------------------------------
# 加载必要的包
library(boot)

# 故障间隔时间数据
x <- aircondit[1]

# 定义统计量函数
meant <- function(x, i) {
  # 计算 bootstrap 样本的平均值
  return(mean(as.matrix(x[i,])))
}

# 执行 bootstrap 过程
set.seed(123) # 设置随机种子以确保结果可重复
b <- boot(x, statistic = meant, R = 2000)

# 计算置信区间
ci_results <- boot.ci(b, type = c("norm", "perc", "basic", "bca"))

# 输出结果
print(ci_results)

##  比较这些区间并解释它们可能不同的原因

##  1. 正态区间与百分位数区间的差异：由于重复样本的分布不是近似正态的，因此正态区间和百分位数区间会有所不同。从重复样本的直方图来看，重复样本的分布是偏斜的。尽管我们正在估计一个平均值，但样本量太小，无法通过中心极限定理（CLT）在这里给出一个好的近似。
##  2. BCa区间的特点：BCa区间是一种百分位数类型的区间，但它同时调整了偏斜和偏差。这意味着BCa区间在考虑了数据分布的偏斜性之后，提供了一个更为准确的区间估计。
##  3. 样本量的影响：样本量对于估计的准确性有重要影响。在本例中，样本量可能不足以使CLT适用，因此正态近似可能不准确。这可能导致正态区间和百分位数区间之间的差异。
##  4. 偏斜性的影响：数据的偏斜性也会影响区间估计。百分位数区间不依赖于分布的正态性，因此它可能在偏斜分布中提供更稳定的估计。而BCa区间则进一步考虑了偏斜性和偏差，因此在这种情况下可能提供最可靠的估计。
##  总结来说，不同的区间估计方法对于数据分布的假设不同，因此在面对偏斜分布和样本量较小的情况时，它们可能会给出不同的结果。BCa区间通过调整偏斜性和偏差，通常在这种情况下提供更为稳健的估计。

# 绘制直方图
hist(b$t, prob = TRUE, main = "")
original_mean <- meant(x, 1:nrow(x))
points(original_mean, 0, cex = 2, pch = 16)

