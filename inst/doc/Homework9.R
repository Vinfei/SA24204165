## -----------------------------------------------------------------------------
# 定义R版本的采样器
gibbsSamplerR <- function(a, b, n, N, x_init) {
  # 初始化
  chain <- matrix(0, nrow = N, ncol = 2)
  x <- x_init
  
  for (i in 1:N) {
    # 从 Beta(x + a, n - x + b) 分布采样 y
    y <- rbeta(1, x + a, n - x + b)
    
    # 从 Binomial(n, y) 分布采样 x
    x <- rbinom(1, n, y)
    
    # 存储结果
    chain[i, ] <- c(x, y)
  }
  
  return(chain)
}

library(Rcpp)

# Rcpp 实现的 Gibbs 采样器
cppFunction('
NumericMatrix gibbsSampler(int a, int b, int n, int N, int x_init) {
  NumericMatrix chain(N, 2);
  int x = x_init;
  double y;
  
  for (int i = 0; i < N; i++) {
    y = R::rbeta(x + a, n - x + b);
    x = R::rbinom(n, y);
    chain(i, 0) = x;
    chain(i, 1) = y;
  }
  
  return chain;
}
')

# 参数设置
a <- 2
b <- 3
n <- 10
N <- 1000
x_init <- 5

# 使用 R 和 Rcpp 实现采样
samples_R <- gibbsSamplerR(a, b, n, N, x_init)
samples_Rcpp <- gibbsSampler(a, b, n, N, x_init)

# 提取 x 和 y 的采样结果
x_R <- samples_R[, 1]
y_R <- samples_R[, 2]
x_Rcpp <- samples_Rcpp[, 1]
y_Rcpp <- samples_Rcpp[, 2]

# (3) 使用 qqplot 比较生成的随机数
par(mfrow = c(1, 2))  # 设置图形布局为 1 行 2 列

# 比较 x 的采样
qqplot(x_R, x_Rcpp, main = "QQ Plot of x", xlab = "R Samples (x)", ylab = "Rcpp Samples (x)")
abline(0, 1, col = "red")

# 比较 y 的采样
qqplot(y_R, y_Rcpp, main = "QQ Plot of y", xlab = "R Samples (y)", ylab = "Rcpp Samples (y)")
abline(0, 1, col = "red")


## -----------------------------------------------------------------------------
library(microbenchmark)

# 定义R版本的采样器
gibbsSamplerR <- function(a, b, n, N, x_init) {
  chain <- matrix(0, nrow = N, ncol = 2)
  x <- x_init
  for (i in 1:N) {
    y <- rbeta(1, x + a, n - x + b)
    x <- rbinom(1, n, y)
    chain[i, ] <- c(x, y)
  }
  return(chain)
}

# Rcpp版本的Gibbs采样器
cppFunction('
NumericMatrix gibbsSampler(int a, int b, int n, int N, int x_init) {
  NumericMatrix chain(N, 2);
  int x = x_init;
  double y;
  
  for (int i = 0; i < N; i++) {
    y = R::rbeta(x + a, n - x + b);
    x = R::rbinom(n, y);
    chain(i, 0) = x;
    chain(i, 1) = y;
  }
  
  return chain;
}
')

# 参数设置
a <- 2
b <- 3
n <- 10
N <- 10000  # 设置较大的采样次数以观察性能差异
x_init <- 5

# 比较两种实现的运行时间
benchmark_results <- microbenchmark(
  R = gibbsSamplerR(a, b, n, N, x_init),
  Rcpp = gibbsSampler(a, b, n, N, x_init),
  times = 10  # 重复测试 10 次
)

# 打印结果
print(benchmark_results)

# 可视化结果
boxplot(benchmark_results, main = "Computation Time Comparison", ylab = "Time (ms)")


