library(gbeta)
library(Runuran)

c = 2; d = 3; kappa = 4; tau = 5; scale = 1
unuran_gbetap_pinv <- pinv.new(
  cdf = function(x) pgbetap(x, c, d, kappa, tau, scale),
  lb = 0, ub = Inf
)

c = 2; d = 3; kappa = 4; tau = 5; 

unuran_gbeta_pinv <- function(){
  pinv.new(
    cdf = function(x) pgbeta(x, c, d, kappa, tau),
    lb = 0, ub = 1
  )
}

dd <- uq(unuran_gbeta_pinv(), 0.5)

