library(gbp)

c <- 2; d <- 3; kappa <- 0.5#-c-d+0.5; sampler pour kappa>0 uniquement
tau <- 1/5

# curve(dgbp(x, c, d, kappa, tau), to=4, col = "red", 
#       lty = "dashed", lwd = 2)

nsims <- 200000
sims <- rgbp(nsims, c, d, kappa, tau)
plot(density(sims, from = 0, to = quantile(sims, 0.95)))
curve(dgbp(x, c, d, kappa, tau), add = TRUE, col = "red", 
      lty = "dashed", lwd = 2)

q <- seq(0.1, 4, length.out = 10)
cdfValues <- pgbp(q, c, d, kappa, tau)
empirical_cdf <- ecdf(sims)
plot(empirical_cdf, xlim = c(0,4))
points(q, cdfValues, pch=19)