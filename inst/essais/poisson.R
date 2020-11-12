T = 15; S = 10
a = 2; b = 3
rho <- S/(T+b)
tau = 4
xo = 5; yo = 12

upost <- function(phi){
  phi^(c+xo-1)*(1+rho*phi)^(-a-xo-yo) / (1+phi/tau)^(c+d)
}

C <- integrate(upost, 0, 4)$value

post <- function(phi) upost(phi) / C

curve(post, from = 0, to = 3)
curve(dgbetap(x, c+xo, a+d+yo, c+d, 1/(rho*tau), scale = 1/rho), add = TRUE, col = "red")

