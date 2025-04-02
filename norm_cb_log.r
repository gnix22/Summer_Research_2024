# defining the density function of our defined ncbl function
pdf_ncbl <- function(x, mu, sig, alpha, lamb) {
f_r <- (2*atanh(1-2*lamb)*(lamb^x)*((1-lamb)^(1-x)))/(1-(2*lamb)) 
F_r <- ((lamb^x) * ((1-lamb)^(1-x))+lamb-1)/((2*lamb)-1)
s_r <- 1 - F_r
z <- exp(-(alpha*log(F_r/(1-F_r))-mu)^2/(2*sig^2))  
numerator <- z * alpha * (s_r + F_r) * f_r # Split numerator and denominator 
denom <- sqrt(2*pi) * s_r * F_r * sig      # for ease of error checking
expr <- numerator / denom
}

# defining distribution function
cdf_ncbl <- function(x, mu, sig, alpha, lamb) {
  cdf_cb <- (lamb ^ x * ((1 - lamb) ^ (1 - x)) + lamb - 1) / (2 * lamb - 1)
  quant_log <- alpha * log(cdf_cb / (1 - cdf_cb))
  cdf_norm <- (1/2) * (1 + erf(quant_log - mu / (sig * sqrt(2))))
}

x1 <- seq(0,1, by= 0.00001) # 0.00001 for smoothing

plot(x1, pdf_ncbl(x1, 0.3, 0.4, 0.9, 0.4),ylim=c(0,5),xlim = c(0, 1),
    ylab = "f(x)", xlab = "", type = "l", lty = 1, col = "black", xaxs="i",yaxs="i")

# let there be lines 
lines(x1, pdf_ncbl(x1, 0.1, 0.9, 0.9, 0.1), type = "l", lty = 2, col = "#E1981A")
lines(x1, pdf_ncbl(x1, 0.1, 0.9, 0.9, 0.9), type = "l", lty = 3, col = "#7e4c0c")
lines(x1, pdf_ncbl(x1, 0.1, 0.9, 0.474, 0.68), type = "l", lty = 4, col = "#e64821")
lines(x1, pdf_ncbl(x1, 0.1, 0.9, 0.5, 0.1), type = "l", lty = 5, col = "#662e07")
lines(x1, pdf_ncbl(x1, 0.1, 0.5, 0.9, 0.1), type = "l", lty = 6, col= "#93360B")
lines(x1, pdf_ncbl(x1, 0.1, 0.5, 0.9, 0.9), type = "l", lty = 7, col="#7D1E08")

col <- c("black", "#E1981A", "#7e4c0c", "#e64821", "#662e07","#93360B","#7D1E08")

legend(0.05,4.9,lwd=2,c(
  expression(mu==0.3~~~~sigma==0.4~~~~alpha==0.9~~~~lambda==0.4),
  expression(mu==0.1~~~~sigma==0.9~~~~alpha==0.9~~~~lambda==0.1),
  expression(mu==0.1~~~~sigma==0.9~~~~alpha==0.9~~~~lambda==0.9),
  expression(mu==0.1~~~~sigma==0.9~~~~alpha==0.474~~~~lambda==0.68),
  expression(mu==0.1~~~~sigma==0.9~~~~alpha==0.5~~~~lambda==0.1),
  expression(mu==0.1~~~~sigma==0.5~~~~alpha==0.9~~~~lambda==0.1),
  expression(mu==0.1~~~~sigma==0.5~~~~alpha==0.9~~~~lambda==0.9)),
 
  lty=1:7,cex=.65, col=col,bty="n")
mtext(text="x", side=1, line=1.75)

