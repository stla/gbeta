library(gbp)
library(Runuran)

ungbpDetails <- function(c, d, kappa, tau){
  distr <- unuran.cont.new(
    cdf = function(q) pgbp(q, c, d, kappa, tau),
    pdf = function(phi) dgbp(phi, c, d, kappa, tau, log = TRUE),
    dpdf = function(x) (c-1)/x + (kappa-c-d)/(1+x) - kappa*tau/(1+tau*x),
    islog = TRUE,
    lb = 0, ub = Inf
  )
  unuran <- unuran.new(distr, method = "arou")
  unuran.details(unuran) 
}


rgbp2 <- function(n, c, d, kappa, tau){
  distr <- unuran.cont.new(
    cdf = function(q) pgbp(q, c, d, kappa, tau),
    pdf = function(phi) dgbp(phi, c, d, kappa, tau, log = TRUE),
    dpdf = function(x) (c-1)/x + (kappa-c-d)/(1+x) - kappa*tau/(1+tau*x),
    islog = TRUE,
    lb = 0, ub = Inf
  )
  unuran <- unuran.new(distr, method = "arou")
  ur(unuran, n) # on peut faire uq ? 
}

c <- 2; d <- 3; kappa <- 0.5; tau <- 1/5 # ok for any kappa
nsims <- 200000
sims <- rgbp2(nsims, c, d, kappa, tau)
plot(density(sims, from = 0, to = quantile(sims, 0.95)))
curve(dgbp(x, c, d, kappa, tau), add = TRUE, col = "red", 
      lty = "dashed", lwd = 2)


library(microbenchmark)
microbenchmark(
  mixture = rgbp(500000, c, d, kappa, tau),
  arou    = rgbp2(500000, c, d, kappa, tau),
  times = 15, 
  unit = "relative"
)
