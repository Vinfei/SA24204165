## -----------------------------------------------------------------------------
# 设置随机种子
set.seed(123)

# 定义标准柯西分布的密度函数
cauchy_pdf <- function(x) {
  1 / (pi * (1 + x^2))
}

# Metropolis-Hastings采样器
mh_sampler <- function(num_samples, proposal_sd) {
  # 初始化参数
  current_value <- rnorm(1, 0, proposal_sd)  # 初始值，从正态分布中采样
  sample_chain <- numeric(num_samples)  # 存储采样结果的向量
  
  # 采样过程
  for (i in 2:num_samples) {
    # 提议分布（正态分布）
    proposed_value <- rnorm(1, current_value, proposal_sd)
    
    # 计算接受概率
    acceptance_ratio <- cauchy_pdf(proposed_value) * dnorm(current_value, proposed_value, proposal_sd) /
                         (cauchy_pdf(current_value) * dnorm(proposed_value, current_value, proposal_sd))
    
    # 接受或拒绝提议
    uniform_random <- runif(1)
    if (uniform_random <= acceptance_ratio) {
      current_value <- proposed_value
    }
    
    # 存储结果
    sample_chain[i] <- current_value
  }
  
  return(sample_chain)
}

# 参数设置
num_samples <- 10000
proposal_sd <- 3
burn_in <- 1000

# 生成样本
sample_chain <- mh_sampler(num_samples, proposal_sd)

# 计算分位数
decile_probs <- seq(0.1, 0.9, 0.1)
post_burn_in_samples <- sample_chain[(burn_in + 1):num_samples]
sample_deciles <- quantile(post_burn_in_samples, decile_probs)
cauchy_deciles <- qcauchy(decile_probs)

# 打印分位数
print(round(rbind(sample_deciles, cauchy_deciles), 3))

# 计算上尾分位数
upper_tail_probs <- seq(0.95, 1, 0.01)
sample_upper_tail <- quantile(post_burn_in_samples, upper_tail_probs)
cauchy_upper_tail <- qcauchy(upper_tail_probs)

# 打印上尾分位数
print(round(rbind(sample_upper_tail, cauchy_upper_tail), 3))

# 绘制QQ图
all_decile_probs <- ppoints(100)
all_sample_deciles <- quantile(post_burn_in_samples, all_decile_probs)
all_cauchy_deciles <- qcauchy(all_decile_probs)
qqplot(all_cauchy_deciles, all_sample_deciles, cex = 0.5)
abline(0, 1)

## -----------------------------------------------------------------------------
# 加载包
library(coda)

# 设置随机种子
set.seed(123)

# 定义标准柯西分布的密度函数
cauchy_pdf <- function(x) {
  1 / (pi * (1 + x^2))
}

# Metropolis-Hastings采样器
mh_sampler <- function(num_samples, proposal_sd) {
  current_value <- rnorm(1, 0, proposal_sd)  # 初始化参数
  sample_chain <- numeric(num_samples)  # 存储采样结果
  
  for (i in 2:num_samples) {
    proposed_value <- rnorm(1, current_value, proposal_sd)  # 提议分布
    
    # 计算接受概率
    acceptance_ratio <- cauchy_pdf(proposed_value) * dnorm(current_value, proposed_value, proposal_sd) /
                         (cauchy_pdf(current_value) * dnorm(proposed_value, current_value, proposal_sd))
    
    # 接受或拒绝提议
    uniform_random <- runif(1)
    if (uniform_random <= acceptance_ratio) {
      current_value <- proposed_value
    }
    
    sample_chain[i] <- current_value  # 存储结果
  }
  
  return(sample_chain)
}

# 参数设置
num_samples <- 10000
proposal_sd <- 3
burn_in <- 1000
num_chains <- 4  # 链的数量

# 初始化多个链
chains <- list()
for (k in 1:num_chains) {
  chains[[k]] <- mh_sampler(num_samples, proposal_sd)
}

# 计算Gelman-Rubin诊断统计量
mcmc_chains <- lapply(chains, function(chain) mcmc(chain[(burn_in + 1):num_samples]))  # 丢弃burn-in
gelman_results <- gelman.diag(mcmc.list(mcmc_chains), autoburnin = FALSE)

# 打印Gelman-Rubin结果
print(gelman_results)

# 检查收敛性，直到所有R_hat < 1.2
while (max(gelman_results$psrf[,1]) >= 1.2) {
  # 生成新的链并更新
  new_chain <- mh_sampler(num_samples, proposal_sd)
  chains <- append(chains, list(new_chain))
  
  # 更新Gelman-Rubin统计量
  mcmc_chains <- lapply(chains, function(chain) mcmc(chain[(burn_in + 1):num_samples]))
  gelman_results <- gelman.diag(mcmc.list(mcmc_chains), autoburnin = FALSE)
  print(gelman_results)
}

# 计算分位数
post_burn_in_samples <- unlist(lapply(chains, function(chain) chain[(burn_in + 1):num_samples]))
decile_probs <- seq(0.1, 0.9, 0.1)
sample_deciles <- quantile(post_burn_in_samples, decile_probs)
cauchy_deciles <- qcauchy(decile_probs)

# 打印分位数
print(round(rbind(sample_deciles, cauchy_deciles), 3))

# 计算上尾分位数
upper_tail_probs <- seq(0.95, 1, 0.01)
sample_upper_tail <- quantile(post_burn_in_samples, upper_tail_probs)
cauchy_upper_tail <- qcauchy(upper_tail_probs)

# 打印上尾分位数
print(round(rbind(sample_upper_tail, cauchy_upper_tail), 3))

# 绘制QQ图
all_decile_probs <- ppoints(100)
all_sample_deciles <- quantile(post_burn_in_samples, all_decile_probs)
all_cauchy_deciles <- qcauchy(all_decile_probs)
qqplot(all_cauchy_deciles, all_sample_deciles, cex = 0.5)
abline(0, 1)

## -----------------------------------------------------------------------------
# 设置参数
N <- 10000  # 总迭代次数
burn <- 2000 
a <- 3 
b <- 5 
n <- 10  # 总的试验次数

# 初始化
x <- numeric(N)
y <- numeric(N)
x[1] <- rbinom(1, n, 0.5)  # 随机初始化x
y[1] <- rbeta(1, x[1] + a, n - x[1] + b)  # 随机初始化y

# 吉布斯采样
for (i in 2:N) {
  x[i] <- rbinom(1, n, y[i - 1])  # 从条件分布中采样x
  y[i] <- rbeta(1, x[i] + a, n - x[i] + b)  # 从条件分布中采样y
}
xb <- x[(burn + 1):N]

# 计算边际分布的估计值
f1 <- table(xb) / length(xb)

# 计算真实的边际概率质量函数
i <- 0:n
fx <- choose(n, i) * beta(i + a, n - i + b) / beta(a, b)

# 绘制结果
barplot(fx, space = 0, ylim = c(0, 0.15), xlab = "n", main = "p(n)=bar; est=points")
points(0:n + 0.5, f1)

## -----------------------------------------------------------------------------
# 加载包
library(coda)

# 设置参数
N <- 10000    
burn <- 2000  
a <- 3
b <- 5
n <- 10     
num_chains <- 4  # 链的数量

# 初始化多个链
chains <- list()
for (k in 1:num_chains) {
  x <- numeric(N)
  y <- numeric(N)
  x[1] <- rbinom(1, n, 0.5)  # 每条链独立随机初始化x
  y[1] <- rbeta(1, x[1] + a, n - x[1] + b)  # 随机初始化y
  
  # 吉布斯采样
  for (i in 2:N) {
    x[i] <- rbinom(1, n, y[i - 1])
    y[i] <- rbeta(1, x[i] + a, n - x[i] + b)
  }
  
  # 将结果存入链列表
  chains[[k]] <- x[(burn + 1):N]
}

# 计算Gelman-Rubin诊断统计量
# 将每条链转换为mcmc对象
mcmc_chains <- lapply(chains, mcmc)

# 使用gelman.diag函数计算R_hat值
gelman_results <- gelman.diag(mcmc.list(mcmc_chains), autoburnin = FALSE)

# 显示R_hat值
print(gelman_results)
# 重复上面的采样过程并计算R_hat，直到所有R_hat < 1.2

# 可视化边际分布的估计值
# 取所有链的样本合并，计算边际分布估计值
all_samples <- unlist(chains)
f1 <- table(all_samples) / length(all_samples)

# 计算真实的边际概率质量函数
i <- 0:n
fx <- choose(n, i) * beta(i + a, n - i + b) / beta(a, b)

# 绘制结果
barplot(fx, space = 0, ylim = c(0, 0.15), xlab = "n", main = "p(n)=bar; est=points")
points(0:n + 0.5, f1)


