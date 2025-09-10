library(stats)
library(stats4)
library(numDeriv)
library(pracma)

##### ---- ORIGINAL CODE FOR REFERENCE ---- #####
#set.seed(15) # seed of original testing
# load in data (original code)
#df <- read.csv("nearest-earth-objects(1910-2024).csv")
#diameter_col <- df[, c(5)]
# based on five no summary and fences, upper inner bound for data to be within support.
#dataset <- diameter_col[!is.na(diameter_col) & diameter_col<0.7126185 & diameter_col>0] # & diameter_col<1 & diameter_col>0.1
#data <- sample(dataset, 200)

# load rds files
data <- readRDS("data/sampled_neo_4_4_2.rds")

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
distributions <- list(
  # ncbl list of elements
  ncbl = list(
    name = "NCBL",
    pdf = function(params, x) {
      lambda <- params[1]; alpha <- params[2]; mu <- params[3]; sigma <- params[4]
      f_r <- (2 * atanh(1 - 2 * lambda) * (lambda^x) * ((1 - lambda)^(1 - x))) / (1 - 2 * lambda)
      F_r <- ((lambda^x) * ((1 - lambda)^(1 - x)) + lambda - 1) / ((2 * lambda) - 1)
      s_r <- 1 - F_r
      z <- exp(-(alpha * log(F_r / s_r) - mu)^2 / (2 * sigma^2))
      num <- z * alpha * (s_r + F_r) * f_r
      denom <- sqrt(2 * pi) * s_r * F_r * sigma
      num / denom
    },
    cdf = function(params, x) {
      lambda <- params[1]; alpha <- params[2]; mu <- params[3]; sigma <- params[4]
      F_r <- ((lambda^x) * ((1 - lambda)^(1 - x)) + lambda - 1) / ((2 * lambda) - 1)
      s_r <- 1 - F_r
      0.5 * (1 + erf((alpha * log(F_r / s_r) - mu) / (sigma * sqrt(2))))
    },
    log_likelihood = function(params, x) {
      f <- distributions$ncbl$pdf(params, x)
      -sum(log(f))
    },
    init = c(0.1, 0.5, 0.2, 0.9),
    lower = c(1e-10, 1e-10, -Inf, 1e-10),
    upper = c(1 - 1e-10, Inf, Inf, Inf),
    method = "L-BFGS-B"
  ),
  # cb list of elements
  cb = list(
    name = "CB",
    pdf = function(params, x) {
      lambda <- params[1]
      (2 * atanh(1 - 2 * lambda) * (lambda^x) * ((1 - lambda)^(1 - x))) / (1 - 2 * lambda)
    },
    cdf = function(params, x) {
      lambda <- params[1]
      ((lambda^x) * ((1 - lambda)^(1 - x)) + lambda - 1) / ((2 * lambda) - 1)
    },
    log_likelihood = function(params, x) {
      f <- distributions$cb$pdf(params, x)
      -sum(log(f))
    },
    init = c(0.05),
    lower = c(0.01),
    upper = c(0.99),
    method = "Brent"
  ),
  exp = list(
    name = "EXP",
    pdf = function(params, x) {
      alpha <- params[1]
      alpha * exp(-alpha * x)
    },
    cdf = function(params, x) {
      alpha <- params[1]
      1 - exp(-alpha * x)
    },
    log_likelihood = function(params, x) {
      f <- distributions$exp$pdf(params, x)
      -sum(log(f))
    },
    init = c(0.1),
    lower = c(0.1),
    upper = c(Inf),
    method = "L-BFGS-B"
  ),
  # beta list of elements
  beta = list(
    name = "BETA",
    pdf = function(params, x) {
      alpha <- params[1]; beta <- params[2]
      dbeta(x, alpha, beta)
    },
    cdf = function(params, x) {
      alpha <- params[1]; beta <- params[2]
      pbeta(x, alpha, beta)
    },
    log_likelihood = function(params, x) {
      f <- distributions$beta$pdf(params, x)
      -sum(log(f))
    },
    init = c(1.9, 4.5),
    lower = c(1e-5, 1e-5),
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

  result <- optim(
    par = dist$init,
    fn = function(p) dist$log_likelihood(p, data),
    gr = grad_fn,
    method = dist$method,
    lower = dist$lower,
    upper = dist$upper,
    control = list(maxit = 1000, trace = 3)
  )

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
pdf("neo_plot.pdf", width = 10)
hist(data, freq = FALSE, main = "Estimated Minimum Diameters",
     xlab = "Distance (in km)", ylab = "Density", xlim = c(0, 1), ylim = c(0, 4.75))
x1 <- seq(0, 1, by = 1e-4)
colors <- c("#000000", "#566E3D", "#FF934F", "#D30C7B")
lt <- c(1, 4, 5, 6)
i <- 1
for (name in names(distributions)) {
  lines(x1, distributions[[name]]$pdf(results[[name]]$par, x1), col = colors[i], lty = lt[i], lwd = 2)
  i <- i + 1
}
legend(0.75, 2, legend = toupper(names(distributions)), col = colors, lty = lt, lwd = 2)
dev.off()