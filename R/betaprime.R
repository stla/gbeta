#' @importFrom stats rbeta dbeta pbeta
NULL

rbetaprime <- function(n, c, d, lambda){
  stopifnot(lambda > 0)
  u <- rbeta(n, c, d)
  lambda * u/(1 - u)
}

dbetaprime <- function(x, c, d, lambda, log){
  stopifnot(lambda > 0)
  if(log){ # TODO: handle x < 0 - DONE in dgbetap
    log(lambda) - 2*log(lambda+x) + dbeta(x/(lambda + x), c, d, log = TRUE)
  }else{
    lambda/(lambda + x)^2 * dbeta(x/(lambda + x), c, d)
  }
}

pbetaprime <- function(q, c, d, lambda){
  stopifnot(lambda > 0)
  pbeta(q/(lambda + q), c, d)
}