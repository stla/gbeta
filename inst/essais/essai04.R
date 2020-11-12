library(gbp)

c <- 2; d <- 3; kappa <- 0.5; tau <- 1/5
nsims <- 200000
sims_mixture <- rgbeta(nsims, c, d, kappa, tau, method = "mixture")
sims_arou <- rgbeta(nsims, c, d, kappa, tau, method = "arou")
plot(density(sims_mixture, from = 0, to = 1), xlim = c(0, 1))
lines(density(sims_arou, from = 0, to = 1), col = "blue", lty = "dashed")
curve(dgbeta(x, c, d, kappa, tau), add = TRUE, col = "red", 
      lty = "dashed", lwd = 2)

plot(ecdf(sims_mixture), xlim = c(0, 1))
lines(ecdf(sims_arou), col = "blue", lty = "dashed")
curve(pgbeta(x, c, d, kappa, tau), add = TRUE, col = "red", 
      lty = "dashed", lwd = 2)

library(microbenchmark)
microbenchmark(
  mixture = rgbp(500000, c, d, kappa, tau, "mixture"),
  arou    = rgbp(500000, c, d, kappa, tau, "arou"),
  times = 15, 
  unit = "relative"
)
