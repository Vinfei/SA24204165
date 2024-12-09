#' @import knitr
#' @import ggplot2
#' @import boot
#' @import bootstrap
#' @import DAAG
#' @import MASS
#' @import coda
#' @import Rcpp
#' @import microbenchmark
#' @import truncnorm
#' @import stats
#' @import glmnet
#' @import expm
#' @import flare

#' @title High-dimensional Regression with LASSO and Confidence Intervals
#' @description Implements a high-dimensional regression procedure using LASSO, computes confidence intervals, and provides unbiased estimates.
#' @param X A numeric matrix of predictors.
#' @param y A numeric vector of response values.
#' @param alpha Significance level for confidence intervals (default: 0.05).
#' @param lambda Regularization parameter for LASSO (default: NULL, calculated internally).
#' @param max_iter Maximum number of iterations for optimization algorithms (default: 50).
#' @param tol Convergence tolerance for optimization (default: 1e-2).
#' @param verbose Logical; if TRUE, displays progress (default: TRUE).
#' @return A list containing:
#'   \item{estimates}{LASSO regression coefficients.}
#'   \item{unbiased_estimates}{Unbiased regression coefficients.}
#'   \item{confidence_intervals}{Matrix of confidence intervals for the coefficients.}
#'   \item{noise_sd}{Estimated standard deviation of noise.}
#' @examples
#' {
#' set.seed(123)
#' n <- 100
#' p <- 50
#' X <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(rep(2, 5), rep(0, p - 5))
#' y <- X %*% beta_true + rnorm(n)
#' result <- SSLASSO(X, y, alpha = 0.05)
#' print(result$estimates)
#' print(result$unbiased_estimates)
#' print(result$confidence_intervals)
#' }
#' @name SSLASSO
#' @export
SSLASSO <- function(X, y, alpha = 0.05, lambda = NULL, max_iter = 50, tol = 1e-2, verbose = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  scaled_X <- scale(X, scale = FALSE)
  design_norm <- sqrt(colSums(scaled_X^2) / n)
  scaled_X <- sweep(scaled_X, 2, design_norm, `/`)
  
  # Get Lasso regression coefficients
  lasso_fit <- Lasso(scaled_X, y, lambda = lambda)
  
  # Ensure lasso_fit is a vector with length equal to p
  if (length(lasso_fit) != p) {
    stop("Lasso() output should be a vector of length p.")
  }
  
  cov_matrix <- crossprod(scaled_X) / n
  inverse_matrix <- compute_inverse(cov_matrix, n, max_iter = max_iter, tol = tol, verbose = verbose)
  unbiased_estimates <- lasso_fit + inverse_matrix %*% crossprod(scaled_X, y - scaled_X %*% lasso_fit) / n
  
  weight_matrix <- inverse_matrix %*% cov_matrix %*% inverse_matrix
  noise_est <- Noise(unbiased_estimates, weight_matrix, n)
  
  interval_half_width <- qnorm(1 - alpha / 2) * noise_est$noise_sd * sqrt(diag(weight_matrix)) / sqrt(n)
  
  list(
    estimates = lasso_fit,
    unbiased_estimates = unbiased_estimates,
    confidence_intervals = cbind(
      unbiased_estimates - interval_half_width,
      unbiased_estimates + interval_half_width
    ),
    noise_sd = noise_est$noise_sd
  )
}

#' @title Soft Threshold Function
#' @description Applies a soft-thresholding operation to a given value.
#' @param value Numeric; the input value to threshold.
#' @param threshold Numeric; the threshold value.
#' @return A numeric value after applying the soft-threshold.
#' @examples
#' {
#' result <- soft_threshold(3.5, 2)
#' print(result)
#' }
#' @name soft_threshold
#' @export
soft_threshold <- function(value, threshold) {
  if (value > threshold) {
    return(value - threshold)
  } else if (value < -threshold) {
    return(value + threshold)
  } else {
    return(0)
  }
}

#' @title Solve One Row for Constrained Optimization
#' @description Solves a single row of a covariance matrix for constrained optimization in high-dimensional regression.
#' @param cov_matrix Numeric matrix; the covariance matrix.
#' @param row_idx Integer; index of the row to solve.
#' @param constraint Numeric; the constraint parameter.
#' @param max_iter Integer; maximum number of iterations (default: 50).
#' @param tol Numeric; convergence tolerance (default: 1e-2).
#' @return A list containing:
#'   \item{solution}{The solution vector for the row.}
#'   \item{iterations}{Number of iterations performed.}
#' @examples
#' {
#' cov_matrix <- matrix(c(4, 2, 2, 3), nrow = 2)
#' result <- solve_one_row(cov_matrix, 1, 0.1)
#' print(result)
#' }
#' @name solve_one_row
#' @export
solve_one_row <- function(cov_matrix, row_idx, constraint, max_iter = 50, tol = 1e-2) {
  n <- nrow(cov_matrix)
  rho <- max(abs(cov_matrix[row_idx, -row_idx])) / cov_matrix[row_idx, row_idx]
  initial_constraint <- rho / (1 + rho)
  coeff <- numeric(n)
  
  if (constraint >= initial_constraint) {
    coeff[row_idx] <- (1 - initial_constraint) / cov_matrix[row_idx, row_idx]
    return(list(solution = coeff, iterations = 0))
  }
  
  residual <- -cov_matrix %*% coeff
  residual[row_idx] <- residual[row_idx] + 1
  iter <- 1
  prev_coeff <- coeff
  
  while (iter <= max_iter) {
    for (j in 1:n) {
      residual[j] <- residual[j] + coeff[j] * cov_matrix[j, j]
      coeff[j] <- soft_threshold(residual[j], constraint) / cov_matrix[j, j]
      residual[j] <- residual[j] - coeff[j] * cov_matrix[j, j]
    }
    
    if (sqrt(sum((coeff - prev_coeff)^2)) < tol * sqrt(sum(coeff^2))) {
      break
    }
    
    prev_coeff <- coeff
    iter <- iter + 1
  }
  
  list(solution = coeff, iterations = iter)
}

#' @title Compute Sparse Inverse Matrix
#' @description Computes the sparse inverse of a covariance matrix using row-wise constrained optimization.
#' @param cov_matrix Numeric matrix; the covariance matrix.
#' @param sample_size Integer; the number of samples.
#' @param constraint Numeric; constraint parameter (default: NULL, calculated internally).
#' @param max_iter Integer; maximum number of iterations (default: 50).
#' @param tol Numeric; convergence tolerance (default: 1e-2).
#' @param verbose Logical; if TRUE, displays progress (default: TRUE).
#' @return A numeric matrix representing the sparse inverse of the input matrix.
#' @examples
#' {
#' cov_matrix <- matrix(c(4, 2, 2, 3), nrow = 2)
#' result <- compute_inverse(cov_matrix, sample_size = 100)
#' print(result)
#' }
#' @name compute_inverse
#' @export
compute_inverse <- function(cov_matrix, sample_size, constraint = NULL, max_iter = 50, tol = 1e-2, verbose = TRUE) {
  p <- nrow(cov_matrix)
  inverse_matrix <- matrix(0, nrow = p, ncol = p)
  
  if (is.null(constraint)) {
    constraint <- qnorm(1 - 0.1 / p^2) / sqrt(sample_size)
  }
  
  for (i in 1:p) {
    if (verbose && i %% (p / 10) == 0) {
      cat(sprintf("%d%% completed\n", round(i / p * 100)))
    }
    
    row_solution <- solve_one_row(cov_matrix, i, constraint, max_iter, tol)
    inverse_matrix[i, ] <- row_solution$solution
  }
  
  inverse_matrix
}

#' @title Fit LASSO Model
#' @description Fits a LASSO regression model using either glmnet or square-root Lasso.
#' @param X A numeric matrix of predictors.
#' @param y A numeric vector of response values.
#' @param lambda Regularization parameter for LASSO (default: NULL).
#' @param intercept Logical; whether to include an intercept in the model (default: TRUE).
#' @return A numeric vector of LASSO regression coefficients.
#' @examples
#' {
#' X <- matrix(rnorm(100 * 10), nrow = 100)
#' y <- rnorm(100)
#' result <- Lasso(X, y)
#' print(result)
#' }
#' @name Lasso
#' @export
Lasso <- function(X, y, lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)

  if (is.null(lambda)) {
    lambda <- sqrt(qnorm(1 - (0.1 / p)) / n)
  }

  outLas <- glmnet(X, y, family = "gaussian", alpha = 1, intercept = intercept, lambda = lambda)
  
  if (intercept) {
    coef_out <- as.vector(coef(outLas, s = lambda))
    return(coef_out[-1])  
  } else {
    return(as.vector(coef(outLas, s = lambda))[-1])  
  }
}

#' @title Estimate Noise Standard Deviation
#' @description Estimates the noise standard deviation using the residuals.
#' @param yh Numeric vector; the predicted values or regression coefficients.
#' @param A Numeric matrix; the covariance matrix (or similar matrix).
#' @param n Integer; the number of samples.
#' @return A list containing:
#'   \item{noise_sd}{Estimated noise standard deviation.}
#'   \item{nz}{Number of non-zero coefficients or significant coefficients.}
#' @name Noise
#' @export
Noise <- function(yh, A, n) {
  ynorm <- sqrt(n) * (yh / sqrt(diag(A)))
  sd_hat0 <- mad(ynorm)
  
  zeros <- (abs(ynorm) < 3 * sd_hat0)
  y2norm <- sum(yh[zeros]^2)
  Atrace <- sum(diag(A)[zeros])
  sd_hat1 <- sqrt(n * y2norm / Atrace)
  
  ratio <- sd_hat0 / sd_hat1
  if (max(ratio, 1 / ratio) > 2) {
    print("Warning: Noise estimate problematic")
  }
  
  s0 <- sum(zeros == FALSE)
  return(list("noise_sd" = sd_hat1, "nz" = s0))
}
