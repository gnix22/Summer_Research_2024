library(stats)
library(pracma)
library(stats4)

set.seed(15) # for reproducibility

# define dataframes
df <- read.csv("datasets/nearest-earth-objects(1910-2024).csv")
diameter_col <- df[, c(5)]

df1 <- read.csv("datasets/darwin_data.csv")
df1_main_col <- df1[, c(9)]
df1_col <- df1_main_col[df1_main_col < 0.2903525]

df2 <- read.csv("datasets/abalone_cleaned_col5_data_v2.csv")
df2_main_col <- df2[, c(1)]

dataset <- diameter_col[!is.na(diameter_col) & diameter_col<0.7126185 & diameter_col>0] # & diameter_col<1 & diameter_col>0.1
#sampled <- sample(dataset, 1000, replace = FALSE)

data <- sample(dataset, 200)
#data <- sample(df2_main_col, 100)

# NCBL Dist

#ncbl_pdf, used for modeling optim params
ncbl_pdf <- function(params, x) {
    lambda <- params[1]
    alpha <- params[2]
    mu <- params[3]
    sigma <- params[4]
    cLamb = (2 * atanh(1 - 2 * lambda)) / (1 - 2*lambda)
    f_r = (2 * atanh(1 - 2 * lambda) * (lambda ^ x) *((1 - lambda) ^ (1 - x))) / (1-(2 * lambda)) 
    F_r = ((lambda^x) * ((1 - lambda) ^ (1 - x)) + lambda - 1) / ((2 * lambda) - 1)
    s_r = 1 - F_r
    z = exp(-(alpha * log(F_r /(1 - F_r)) - mu)^2 /(2 * sigma^2))  
    numerator = z * alpha * (s_r + F_r) * f_r
    denom = sqrt(2 * pi) * s_r * F_r * sigma
    f = numerator / denom
}

#ncbl_cdf, used for ks test later
ncbl_cdf <- function(params, x) {
    lambda <- params[1]
    alpha <- params[2]
    mu <- params[3]
    sigma <- params[4]
    f_r <- (2 * atanh(1 - 2 * lambda) * (lambda ^ x) *((1 - lambda) ^ (1 - x))) / (1-(2 * lambda)) 
    F_r <- ((lambda^x) * ((1 - lambda) ^ (1 - x)) + lambda - 1) / ((2 * lambda) - 1)
    s_r <- 1 - F_r
    f <- (1 / 2) * (1 + erf((alpha * log(F_r / s_r) - mu) / (sigma * sqrt(2))))
}

#define log likelihood function
ncbl_log_likelihood <- function(params, x) {
  lambda <- params[1]
  alpha <- params[2]
  mu <- params[3]
  sigma <- params[4]
  
  cLamb = (2 * atanh(1 - 2 * lambda)) / (1 - 2*lambda)
  f_r = (2 * atanh(1 - 2 * lambda) * (lambda ^ x) *((1 - lambda) ^ (1 - x))) / (1-(2 * lambda)) 
  F_r = ((lambda^x) * ((1 - lambda) ^ (1 - x)) + lambda - 1) / ((2 * lambda) - 1)
  s_r = 1 - F_r
  z = exp(-(alpha * log(F_r /(1 - F_r)) - mu)^2 /(2 * sigma^2))  
  numerator = z * alpha * (s_r + F_r) * f_r
  denom = sqrt(2 * pi) * s_r * F_r * sigma
  f = numerator / denom
    
  ll = -sum(log(f))
}

# mle using l-bfgs-b algorithm
bfgs_optimization_ncbl <- function(data) {
  init_params <- c(0.1, 0.5, 0.2, 0.9) 
  
  result <- optim(
    par = init_params,
    fn = ncbl_log_likelihood,
    x = data,
    method = "L-BFGS-B", 
    lower = c(1e-10, 1e-10, -Inf, 1e-10), #lower bounds
    upper = c(1-(1e-10), Inf, Inf, Inf), #upper bounds
    control = list(maxit = 1000, trace = 3),
    hessian = TRUE
  )
  
  return(result)
}

# CB Dist
cb_pdf <- function(params, x) {
    lambda <- params[1]
    cLamb = (2 * atanh(1 - 2 * lambda)) / (1 - 2*lambda)
    f_r = (2 * atanh(1 - 2 * lambda) * (lambda ^ x) *((1 - lambda) ^ (1 - x))) / (1-(2 * lambda))    
}

#cb cdf, used for ks test
cb_cdf <- function(params, x) {
  lambda <- params[1]
  F_r = ((lambda^x) * ((1 - lambda) ^ (1 - x)) + lambda - 1) / ((2 * lambda) - 1)
}

cb_log_likelihood <- function(params, x) {
    lambda <- params[1]
    x = data
    cLamb = (2 * atanh(1 - 2 * lambda)) / (1 - 2*lambda)
    f = (2 * atanh(1 - 2 * lambda) * (lambda ^ x) *((1 - lambda) ^ (1 - x))) / (1-(2 * lambda))
    l = prod(f)
    ll = -log(l)
}  

# MLE
bfgs_optimization_cb <- function(params, x) {
  init_params <- c(0.001) 
  
  result <- optim(
    par = init_params,
    fn = cb_log_likelihood,
    x = data,
    method = "Brent", 
    lower = c(0.01), #lower bounds
    upper = c(0.99), #upper bounds
    control = list(maxit = 1000, trace = 3),
    hessian = TRUE
  )
  
  return(result)
}
# Expo Dist
exp_pdf <- function(params, x) {
  alpha <- params[1]
  f <- alpha * exp(-alpha * x) 
}

exp_cdf <- function(params, x) {
  alpha <- params[1]
  f <- 1 - exp(-alpha * x)
}

exp_log_likelihood <- function(params, x) {
  alpha <- params[1]
  f <- alpha * exp(-alpha * x)
  l <- prod(f)
  ll <- -log(l)
}

bfgs_optimization_exp <- function(params, x) {
  init_params <- c(.1)
  
  result <- optim(
    par = init_params,
    fn = exp_log_likelihood,
    x = data,
    method = "L-BFGS-B",
    lower = c(0.1), #lower bounds 1e-6
    upper = c(Inf), #upper bounds
    control = list(maxit = 1000, trace = 3),
    hessian = TRUE
  )
  
  return(result)
}
# Beta Dist
beta_pdf <- function(params, x) {
  alpha <- params[1]
  beta <- params[2]
  dbeta(x, alpha, beta)
}

beta_pdf_v2 <- function(params, x) {
  alpha <- params[1]
  beta_par <- params[2]
  f <- (1 / beta(alpha, beta_par)) * x ^ (alpha-1) * (1-x) ^ (beta_par - 1)
}

beta_cdf <- function(params, x) {
  alpha <- params[1]
  beta <- params[2]
  pbeta(x, alpha, beta)
}

beta_log_likelihood <- function(params, x) {
  alpha <- params[1]
  beta_par <- params[2]
  f <- (1 / beta(alpha, beta_par)) * x ^ (alpha-1) * (1-x) ^ (beta_par - 1)
  l <- prod(f)
  ll <- -log(l)
}

bfgs_optimization_beta <- function(params, x) {
    init_params <- c(1.9, 4.5)
  
  result <- optim(
    par = init_params,
    fn = beta_log_likelihood,
    x = data,
    method = "L-BFGS-B",
    lower = c(1e-5,1e-5), #lower bounds 1e-5
    upper = c(Inf,Inf), #upper bounds
    control = list(maxit = 1000, trace = 3),
    hessian = TRUE
  )
  
  return(result)
}

##data section
##cleaning data

##run BFGS optimization
optim_res_ncbl <- bfgs_optimization_ncbl(data)
optim_res_cb <- bfgs_optimization_cb(data)
optim_res_exp <- bfgs_optimization_exp(data)
optim_res_beta <- bfgs_optimization_beta(data)

print(optim_res_ncbl)
print(optim_res_cb)
print(optim_res_exp)
print(optim_res_beta)

#get params of optimized params (for modeling)

res_params_ncbl <- optim_res_ncbl$par
res_params_cb <- optim_res_cb$par
res_params_exp <- optim_res_exp$par
res_params_beta <- optim_res_beta$par

epsilon <- 1e-3 # define epsilon for covariance matrix making positive values 

# covariance matrices section
hess_ncbl <- optim_res_ncbl$hessian
epsilon_ncbl <- max(abs(diag(hess_ncbl))) * epsilon  # Set epsilon based on the largest diagonal element
cov_ncbl <- solve(hess_ncbl + diag(epsilon_ncbl, nrow(hess_ncbl)))

cat("COVARIANCE MATRIX OF NCBL\n")
print(cov_ncbl)

cat("STANDARD ERRORS (LAMBDA, ALPHA, MU, SIGMA)\n")
print(diag(sqrt(cov_ncbl)))

hess_cb <- optim_res_cb$hessian
epsilon_cb <- max(abs(diag(hess_cb))) * epsilon  # Set epsilon based on the largest diagonal element
cov_cb <- solve(hess_cb + diag(epsilon_cb, nrow(hess_cb)))
cat("COVARIANCE MATRIX OF CB\n")
print(cov_cb)
cat("STANDARD ERRORS (LAMBDA)\n")
print(diag(sqrt(cov_cb)))

hess_exp <- optim_res_exp$hessian
epsilon_exp <- max(abs(diag(hess_exp))) * epsilon
cov_exp <- solve(hess_exp + diag(epsilon_exp, nrow(hess_exp)))
cat("COVARIANCE MATRIX OF EXP\n")
print(cov_exp)
cat("STANDARD ERRORS (ALPHA)\n")
print(diag(sqrt(cov_exp)))

hess_beta <- optim_res_beta$hessian
epsilon_beta <- max(abs(diag(hess_beta))) * epsilon
cov_beta <- solve(hess_beta + diag(epsilon_beta, nrow(hess_beta)))
cat("COVARIANCE MATRIX OF BETA\n")
print(cov_beta)
cat("STANDARD ERRORS (ALPHA, BETA)\n")
print(diag(sqrt(cov_beta)))
# tests section
# define empirical cdf for use in ks statistic
ecdf_data <- ecdf(data)
x_values <- seq(min(data), max(data), length.out = 1000)
ecdf_values <- ecdf_data(x_values)
n <- length(data)

# for computing p-value of ks test, emulates pkstwo()
pkstwo <- function(x) {
  if (x < 0) return(0)
  return(1 - (1 - x)^(2))
}

# ncbl tests
ncbl_likelihood_val <- exp(optim_res_ncbl$value)
cat("NCBL LIKELIHOOD: \n", ncbl_likelihood_val,"\n")
cat("NCBL AIC\n")

ncbl_aic <- 2 * 4 - 2 * log(ncbl_likelihood_val)
print(ncbl_aic)

cat("NCBL BIC\n")
ncbl_bic <- log(length(data)) * 4 - 2 * log(ncbl_likelihood_val)
print(ncbl_bic)

#ncbl ks test since ks.test won't work properly
ncbl_theoretical_values <- ncbl_cdf(res_params_ncbl, x_values)
ncbl_differences <- abs(ecdf_values - ncbl_theoretical_values)
ncbl_ks <- max(ncbl_differences)
ncbl_pval <- pkstwo(ncbl_ks)
cat("K-S STATISTIC", ncbl_ks, "\n", "P-VAL", ncbl_pval, "\n")

cb_likelihood_val <- exp(optim_res_cb$value)
cat("CB LIKE: ", cb_likelihood_val)
cat("CB AIC\n")
cb_aic <- 2 * 3 - 2 * log(cb_likelihood_val)
print(cb_aic)
cat("CB BIC\n")
cb_bic <- log(length(data)) * 3 - 2 * log(cb_likelihood_val)
print(cb_bic)

# # #ks test for cb 
cb_theoretical_values <- cb_cdf(res_params_cb, x_values)
cb_differences <- abs(ecdf_values - cb_theoretical_values)
cb_ks <- max(cb_differences)
cb_pval <- pkstwo(cb_ks)

cat("K-S STATISTIC", cb_ks, "\n", "P-VAL", cb_pval, "\n")

# Exponential
exp_likelihood_val <- exp(optim_res_exp$value)
cat("EXP LIKE:", exp_likelihood_val)
cat("EXP AIC\n")
exp_aic <- 2 * 3 - 2 * log(exp_likelihood_val)
print(exp_aic)
cat("EXP BIC\n")
exp_bic <- log(length(data)) * 3 - 2 * log(exp_likelihood_val)
print(exp_bic)

exp_theoretical_values <- exp_cdf(res_params_exp, x_values)
exp_differences <- abs(ecdf_values - exp_theoretical_values)
exp_ks <- max(exp_differences)
exp_pval <- pkstwo(exp_ks)

cat("K-S STATISTIC", exp_ks, "\n", "P-VAL", exp_pval, "\n")

# Beta
beta_likelihood_val <- exp(optim_res_beta$value)
cat("BETA LIKE: \n", beta_likelihood_val, "\n")
cat("BETA AIC\n")
beta_aic <- 2 * 3 - 2 * log(beta_likelihood_val)
print(beta_aic)
cat("BETA BIC\n")
beta_bic <- log(length(data)) * 3 - 2 * log(beta_likelihood_val)
print(beta_bic)

# # #ks test for cb 
beta_theoretical_values <- beta_cdf(res_params_beta, x_values)
beta_differences <- abs(ecdf_values - beta_theoretical_values)
beta_ks <- max(beta_differences)
beta_pval <- pkstwo(beta_ks)

cat("K-S STATISTIC", beta_ks, "\n", "P-VAL", beta_pval, "\n")

#modeling section
png(filename = "test.png")
x1 <- seq(0, 1, by=1e-4) # 1e-4 for smoothing
hist(data, freq = FALSE, main = "Estimated Minimum Diameters", 
    xlab = "Distance (in km)", ylab = "Density", xlim = c(0,1), ylim = c(0,4.75)) #, , xaxt = 'n'
lines(x1, ncbl_pdf(res_params_ncbl,x1), col = "blue", type = "l", lty = 4)
lines(x1, beta_pdf(res_params_beta,x1), col = "#ff6f00", type = "l", lty = 1)
lines(x1, cb_pdf(res_params_cb, x1), col = "#44ff00", type = "l", lty = 5)
lines(x1, exp_pdf(res_params_exp, x1), col = "purple", type = "l", lty = 3)
legend(0.75,2, legend = c("NCBL","CB","EXP", "BETA"), 
      col = c("blue","#44ff00", "purple", "#ff6f00"), lty = c(4,5,3,1))
dev.off()
