library(gbp)

c = 2; d = 3; kappa = 4; tau = 5
x = y = 4; n = 10

upost <- function(u){
  dbinom(x, n, u) * dgbeta(u, c, d, kappa, tau)
}

C <- integrate(upost, 0, 1)$value

post <- function(u) upost(u) / C

curve(post)
curve(dgbeta(x, c+y, d+n-y, kappa, tau), add = TRUE, col = "red")
