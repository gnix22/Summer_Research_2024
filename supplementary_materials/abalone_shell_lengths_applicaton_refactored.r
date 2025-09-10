library(stats)
library(stats4)
library(numDeriv)
library(pracma)
##### ---- ORIGINAL CODE FOR REFERENCE ---- #####
#set.seed(4)
#df2 <- read.csv("data/abalone.csv")
#df2_main_col <- df2[, c(2)]
#df2_col <- df2_main_col
#data <- sample(df2_col, 100)

# load rds files
data <- readRDS("data/abalone_shell_lengths_sample_4-4-2.rds")


# x values for plotting and ks testing
x_values <- seq(min(data), max(data), length.out = 1000)
ecdf_data <- ecdf(data)
ecdf_vals <- ecdf_data(x_values)

# pkstwo hard coded because I couldn't get library to work
pkstwo <- function(x) {
  if (x < 0) return(0)
  return(1 - (1 - x)^2)
}

# list of distributions and components
# Distribution definitions
distributions <- list(
  ncbc = list(
    name = "NCBC",
    pdf = function(params, x) {
      lambda <- params[1]; mu <- params[2]; sigma <- params[3]
      F_r <- ((lambda^x) * ((1 - lambda)^(1 - x)) + lambda - 1) / (2 * lambda - 1)
      f_r <- (2 * atanh(1 - 2 * lambda) * (lambda^x) * ((1 - lambda)^(1 - x))) / (1 - 2 * lambda)
      z <- exp(-((tan(pi * (F_r - 0.5)) - mu)^2) / (2 * sigma^2))
      num <- z * pi * sec(pi * (F_r - 0.5)) * f_r
      denom <- sqrt(2 * pi * sigma^2)
      num / denom
    },
    cdf = function(params, x) {
      lambda <- params[1]; mu <- params[2]; sigma <- params[3]
      F_r <- ((lambda^x) * ((1 - lambda)^(1 - x)) + lambda - 1) / (2 * lambda - 1)
      0.5 * (1 + erf((tan(pi * (F_r - 0.5)) - mu) / (sigma * sqrt(2))))
    },
    log_likelihood = function(params, x) {
      lambda <- params[1]; mu <- params[2]; sigma <- params[3]
      F_r <- ((lambda^x) * ((1 - lambda)^(1 - x)) + lambda - 1) / (2 * lambda - 1)
      f_r <- (2 * atanh(1 - 2 * lambda) * (lambda^x) * ((1 - lambda)^(1 - x))) / (1 - 2 * lambda)
      log_lik <- -0.5 * log(2 * pi * sigma^2) - ((tan(pi * (F_r - 0.5)) - mu)^2) / (2 * sigma^2) +
        log(pi * sec(pi * (F_r - 0.5))) + log(f_r)
      -sum(log_lik)
    },
    init = c(0.8, 0.3, 0.2),
    lower = c(0.1, -Inf, 0.1),
    upper = c(0.9, Inf, Inf),
    method = "L-BFGS-B"
  ),

  norm = list(
    name = "Normal",
    pdf = function(params, x) dnorm(x, mean = params[1], sd = params[2]),
    cdf = function(params, x) pnorm(x, mean = params[1], sd = params[2]),
    log_likelihood = function(params, x) -sum(log(dnorm(x, mean = params[1], sd = params[2]))),
    init = c(mean(data), sd(data)),
    method = "BFGS"
  ),

  cauchy = list(
    name = "Cauchy",
    pdf = function(params, x) dcauchy(x, location = params[1], scale = params[2]),
    cdf = function(params, x) pcauchy(x, location = params[1], scale = params[2]),
    log_likelihood = function(params, x) -sum(log(dcauchy(x, location = params[1], scale = params[2]))),
    init = c(median(data), 0.1),
    lower = c(-Inf, 0.001),
    upper = c(Inf, Inf),
    method = "L-BFGS-B"
  )
)

# gradient and hessian helper functions
make_grad_fn <- function(ll_fn, data) {
  function(params) grad(function(p) ll_fn(p, data), params)
}

make_hess_fn <- function(ll_fn, data) {
  function(params) hessian(function(p) ll_fn(p, data), params)
}

# fit distribution using numerical gradient
fit_distribution <- function(dist, data) {
  grad_fn <- make_grad_fn(dist$log_likelihood, data)
  hess_fn <- make_hess_fn(dist$log_likelihood, data)

  optim_args <- list(
    par = dist$init,
    fn = function(p) dist$log_likelihood(p, data),
    gr = grad_fn,
    method = dist$method,
    control = list(maxit = 1000, trace = 3)
  )

  if (dist$method %in% c("L-BFGS-B", "Brent")) {
    optim_args$lower <- dist$lower
    optim_args$upper <- dist$upper
  }

  result <- do.call(optim, optim_args)
  result$hessian <- hess_fn(result$par)
  return(result)
}
# results
evaluate_fit <- function(result, dist_name, data_len, num_params, cdf_func, data, x_vals, ecdf_vals) {
  loglik <- result$value
  aic <- 2 * num_params + 2 * loglik
  bic <- log(data_len) * num_params + 2 * loglik
  epsilon <- max(abs(diag(result$hessian))) * 1e-3
  cov_matrix <- solve(result$hessian + diag(epsilon, nrow(result$hessian)))
  se <- sqrt(diag(cov_matrix))

  theoretical_cdf <- cdf_func(result$par, x_vals)
  ks_stat <- max(abs(ecdf_vals - theoretical_cdf))
  pval <- pkstwo(ks_stat)

  cat(sprintf("=== %s ===\n", dist_name))
  cat("Estimated Parameters:\n"); print(result$par)
  cat("Final Value (NLL):", loglik, "\n")
  cat("Log-Likelihood:", exp(loglik), "\n")
  cat("AIC:", aic, "\n")
  cat("BIC:", bic, "\n")
  cat("K-S Statistic:", ks_stat, "\n")
  cat("P-Value:", pval, "\n")
  cat("Standard Errors:\n"); print(se)
  cat("\n")

  list(aic = aic, bic = bic, ks = ks_stat, pval = pval, se = se, cov = cov_matrix)
}

# testing stuff
results <- lapply(distributions, fit_distribution, data = data)
evaluations <- mapply(
  FUN = evaluate_fit,
  result = results,
  dist_name = names(distributions),
  num_params = sapply(distributions, function(d) length(d$init)),
  cdf_func = lapply(distributions, function(d) d$cdf),
  MoreArgs = list(data_len = length(data), data = data, x_vals = x_values, ecdf_vals = ecdf_vals),
  SIMPLIFY = FALSE
)

# plotting
pdf("shell_lengths_abalone.pdf", width = 10)
hist(data, freq = FALSE, main = "Shell Lengths in Abalone Dataset",
     xlab = "Shell Lengths", ylab = "Density", xlim = c(0, 1), ylim = c(0, 5))
x1 <- seq(0, 1, by = 1e-4)
colors <- c("#000000", "#566E3D", "#FF934F")
lt <- c(1, 4, 6)
i <- 1
for (name in names(distributions)) {
  lines(x1, distributions[[name]]$pdf(results[[name]]$par, x1), col = colors[i], lty = lt[i], lwd = 2)
  i <- i + 1
}
legend(0.15, 4.5, legend = toupper(names(distributions)), col = colors, lty = lt, lwd = 2)
dev.off()