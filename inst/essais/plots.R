library(gbeta)

par(mar = c(2,1,1,1))
layout(cbind(c(1,2),c(3,4)))
curve(dgbeta(x, 4, 12, 1, 2), axes = FALSE, xlab = NA, ylab = NA, lwd = 2)
axis(1)
curve(dgbeta(x, .1, 1, 2, 0.1), axes = FALSE, xlab = NA, ylab = NA, lwd = 2)
axis(1)
curve(dgbeta(x, 4, 12, 10, 0.01), axes = FALSE, xlab = NA, ylab = NA, lwd = 2)
axis(1)
curve(dgbeta(x, 1, 12, 8, 0.01), axes = FALSE, xlab = NA, ylab = NA, lwd = 2)
axis(1)

par(mar = c(2,1,1,1))
layout(cbind(c(1,2),c(3,4)))
curve(dgbetap(x, 4, 12, 1, 2, 5), to = 10, axes = FALSE, xlab = NA, ylab = NA, lwd = 2)
axis(1)
curve(dgbetap(x, .1, 1, 2, 0.1, 5), to = 0.4, axes = FALSE, xlab = NA, ylab = NA, lwd = 2)
axis(1)
curve(dgbetap(x, 4, 12, 10, 0.01), to = 10, axes = FALSE, xlab = NA, ylab = NA, lwd = 2)
axis(1)
curve(dgbetap(x, 1, 12, 8, 0.01, 5), to = 5, axes = FALSE, xlab = NA, ylab = NA, lwd = 2)
axis(1)

curve(pgbetap(x, 4, 12, 10, 0.01), to = 10)

