library(gbp)

c <- 2; d <- 3; kappa <- 0.5; tau <- 1/5; scale <- 2
nsims <- 200000
sims_mixture <- rgbp(nsims, c, d, kappa, tau, scale, method = "mixture")
sims_arou <- rgbp(nsims, c, d, kappa, tau, scale, method = "arou")
plot(density(sims_mixture, from = 0, to = 8), xlim = c(0, quantile(sims_mixture, 0.95)))
lines(density(sims_arou, from = 0, to = 8), col = "blue", lty = "dashed")
curve(dgbp(x, c, d, kappa, tau, scale), add = TRUE, col = "red", 
      lty = "dashed", lwd = 2)

plot(ecdf(sims_mixture), xlim = c(0, quantile(sims_mixture, 0.95)))
lines(ecdf(sims_arou), col = "blue", lty = "dashed")
curve(pgbp(x, c, d, kappa, tau, scale), add = TRUE, col = "red", 
      lty = "dashed", lwd = 2)

library(microbenchmark)
microbenchmark(
  mixture = rgbp(500000, c, d, kappa, tau, "mixture"),
  arou    = rgbp(500000, c, d, kappa, tau, "arou"),
  times = 15, 
  unit = "relative"
)
