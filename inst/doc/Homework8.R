## -----------------------------------------------------------------------------
# 加载包
library(boot)

# 定义约束矩阵 A1
A1 <- rbind(
  c(2, 1, 1),  # 2x + y + z
  c(1, -1, 3)  # x - y + 3z
)

# 定义约束右侧值 b1
b1 <- c(2, 3)

# 定义目标函数的系数向量
a <- c(4, 2, 9)

# 使用 simplex 函数求解最小化问题
result <- simplex(a = a, A1 = A1, b1 = b1, maxi = FALSE)  

# 输出结果
cat("线性规划求解结果：\n")
print(result)

## -----------------------------------------------------------------------------
# 加载 mtcars 数据集
data(mtcars)

# 定义公式列表
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

# 初始化一个空列表存储模型
models_for <- list()

# 使用 for 循环拟合线性模型
for (i in seq_along(formulas)) {
  models_for[[i]] <- lm(formulas[[i]], data = mtcars)
}

# 打印模型摘要
cat("Using for loop:\n")
for (i in seq_along(models_for)) {
  cat("Model", i, ":\n")
  print(summary(models_for[[i]]))
}

# 使用 lapply() 拟合线性模型
models_lapply <- lapply(formulas, function(f) lm(f, data = mtcars))

# 打印模型摘要
cat("\nUsing lapply():\n")
lapply(seq_along(models_lapply), function(i) {
  cat("Model", i, ":\n")
  print(summary(models_lapply[[i]]))
})


## -----------------------------------------------------------------------------
# 创建自举样本
set.seed(123) 

bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

# 使用 for 循环拟合模型
models_for <- list()
for (i in seq_along(bootstraps)) {
  data <- bootstraps[[i]]
  models_for[[i]] <- lm(mpg ~ disp, data = data)
}

# 查看前几个模型
summary(models_for[[1]])

# 定义一个拟合函数
fit_model <- function(data) {
  lm(mpg ~ disp, data = data)
}

# 使用 lapply() 拟合模型
models_lapply <- lapply(bootstraps, fit_model)

# 查看前几个模型
summary(models_lapply[[1]])

## -----------------------------------------------------------------------------
# 创建公式列表
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

# 定义 Rsquared 提取函数
rsq <- function(mod) summary(mod)$r.squared

# 使用 for 循环拟合模型并提取 Rsquared
models_for <- list()
rsq_for <- numeric(length(formulas))
for (i in seq_along(formulas)) {
  models_for[[i]] <- lm(formulas[[i]], data = mtcars)
  rsq_for[i] <- rsq(models_for[[i]])
}
print(rsq_for)

# 使用 lapply() 拟合模型并提取 Rsquared
models_lapply <- lapply(formulas, function(f) lm(f, data = mtcars))
rsq_lapply <- sapply(models_lapply, rsq)
print(rsq_lapply)

## -----------------------------------------------------------------------------
# 定义自举样本
set.seed(123) 

bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), replace = TRUE)
  mtcars[rows, ]
})

# 定义拟合函数
fit_model <- function(data) lm(mpg ~ disp, data = data)

# 使用 for 循环拟合模型并提取 Rsquared
models_boot_for <- list()
rsq_boot_for <- numeric(length(bootstraps))
for (i in seq_along(bootstraps)) {
  models_boot_for[[i]] <- fit_model(bootstraps[[i]])
  rsq_boot_for[i] <- rsq(models_boot_for[[i]])
}
print(rsq_boot_for)

# 使用 lapply() 拟合模型并提取 Rsquared
models_boot_lapply <- lapply(bootstraps, fit_model)
rsq_boot_lapply <- sapply(models_boot_lapply, rsq)
print(rsq_boot_lapply)


## -----------------------------------------------------------------------------
# 模拟非正态数据的 t 检验性能
set.seed(123) 

trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

# 使用匿名函数提取 p 值
p_values_anonymous <- sapply(trials, function(x) x$p.value)
print(head(p_values_anonymous))

# 使用 sapply() 和 [[ 直接提取 p 值
p_values_direct <- sapply(trials, `[[`, "p.value")
print(head(p_values_direct))

## -----------------------------------------------------------------------------
# 创建 lapply() 变体，结合 Map() 和 vapply()
parallel_lapply <- function(FUN, ...) {
  # 获取所有输入
  inputs <- list(...)
  
  # 使用 Map() 并行遍历输入，传递给 FUN
  result <- do.call(Map, c(list(FUN), inputs))  # 显式使用 do.call 来传递多个参数给 FUN
  
  # 使用 vapply() 将结果转化为向量或矩阵
  # 假设每个输入的结果应该是向量，使用 vapply 进行类型保证
  result_vector <- vapply(result, function(x) x, FUN.VALUE = numeric(1))
  
  return(result_vector)
}

# 示例使用
set.seed(123)

# 使用 parallel_lapply() 遍历两个数值向量，计算它们的和
output <- parallel_lapply(function(x, y) x + y, 1:5, 6:10)
print(output)


## -----------------------------------------------------------------------------
# 自定义快速版本的卡方检验函数
fast_chi_sq_test <- function(observed, expected) {
  # 检查输入是否为两个数值向量，并且没有缺失值
  if (length(observed) != length(expected)) {
    stop("The two vectors must have the same length.")
  }
  if (any(is.na(observed)) || any(is.na(expected))) {
    stop("Input vectors should not contain missing values.")
  }
  
  # 计算卡方统计量
  chi_sq_statistic <- sum((observed - expected)^2 / expected)
  
  return(chi_sq_statistic)
}

# 示例：输入两个没有缺失值的数值向量
observed <- c(10, 20, 30, 40)
expected <- c(15, 15, 35, 35)

# 计算卡方统计量
chi_sq_value <- fast_chi_sq_test(observed, expected)
print(chi_sq_value)


## -----------------------------------------------------------------------------
# 自定义快速版本的 table() 函数
fast_table <- function(x, y) {
  # 检查输入是否为两个整数向量，并且没有缺失值
  if (length(x) != length(y)) {
    stop("The two vectors must have the same length.")
  }
  if (any(is.na(x)) || any(is.na(y))) {
    stop("Input vectors should not contain missing values.")
  }
  
  # 计算联合频率分布（交叉表）
  unique_x <- unique(x)
  unique_y <- unique(y)
  
  # 初始化一个零矩阵来存储交叉频率
  result <- matrix(0, nrow = length(unique_x), ncol = length(unique_y),
                   dimnames = list(unique_x, unique_y))
  
  # 填充矩阵
  for (i in seq_along(x)) {
    result[as.character(x[i]), as.character(y[i])] <- result[as.character(x[i]), as.character(y[i])] + 1
  }
  
  return(result)
}

# 自定义快速卡方检验函数
fast_chi_sq_test <- function(observed, expected) {
  # 获取观测频数
  obs_table <- fast_table(observed, expected)
  
  # 计算期望频数
  total <- sum(obs_table)
  row_totals <- rowSums(obs_table)
  col_totals <- colSums(obs_table)
  expected_table <- outer(row_totals, col_totals, FUN = "*") / total
  
  # 计算卡方统计量
  chi_sq_statistic <- sum((obs_table - expected_table)^2 / expected_table)
  
  return(chi_sq_statistic)
}

# 示例：输入两个没有缺失值的整数向量
observed <- c(1, 2, 1, 2, 1, 2, 1, 1, 2)
expected <- c(1, 1, 2, 2, 1, 1, 2, 1, 2)

# 计算卡方统计量
chi_sq_value <- fast_chi_sq_test(observed, expected)
print(chi_sq_value)

