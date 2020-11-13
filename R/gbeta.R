#' Generalized Beta distribution
#' @description Density, distribution function, quantile function, and random 
#'   generation for the generalized Beta distribution.
#'   
#' @param u numeric vector
#' @param q numeric vector of quantiles
#' @param p numeric vector of probabilities
#' @param n positive integer, the desired number of simulations
#' @param c,d,kappa,tau parameters; they must be strictly positive numbers, 
#'   except \code{kappa} which can take any value
#' @param log logical, whether to return the log-density
#' @param method the method of random generation, \code{"mixture"} or 
#'   \code{"arou"}; only a positive \code{kappa} is allowed for the 
#'   \code{"mixture"} method, but this method is faster
#'   
#' @references 
#' \itemize{
#'   \item Marwa Hamza & Pierre Vallois. 
#'     \emph{On Kummer’s distributions of type two and generalized Beta 
#'           distributions}.
#'     Statistics & Probability Letters 118 (2016), pp. 60-69.
#'     <doi:10.1016/j.spl.2016.03.014>
#'   \item James J. Chen & Melvin R. Novick.
#'     \emph{Bayesian Analysis for Binomial Models with Generalized Beta Prior 
#'           Distributions}.
#'     Journal of Educational Statistics 9, No. 2 (1984), pp. 163-175.
#'     <doi:10.3102/10769986009002163>
#' }
#' 
#' @examples library(gbeta)
#' curve(dgbeta(x, 4, 12, 10, 0.01), axes = FALSE, lwd = 2)
#' axis(1)
#' 
#' @importFrom stats qbeta
#' @importFrom Runuran uq
#' 
#' @rdname GBeta
#' @name GBeta
#' @export
dgbeta <- function(u, c, d, kappa, tau, log = FALSE){ 
  stopifnot(c > 0, d > 0, tau > 0)
  out <- numeric(length(u))
  outside01 <- u < 0 | u > 1
  out[outside01] <- ifelse(log, -Inf, 0)
  inside01 <- !outside01
  if(any(inside01)){
    u <- u[inside01]
    out[inside01] <- if(log){
      -2 * log(1-u) + dgbetap(u/(1-u), c, d, kappa, tau, scale = 1, log = TRUE)
    }else{
      1/(1-u)^2 * dgbetap(u/(1-u), c, d, kappa, tau, scale = 1) 
    }
  }
  out
}

#' @rdname GBeta
#' @export
pgbeta <- function(q, c, d, kappa, tau){ 
  stopifnot(c > 0, d > 0, tau > 0)
  out <- numeric(length(q))
  negs <- q <= 0
  sup1 <- q >= 1
  out[sup1] <- 1
  in01 <- !(negs || sup1)
  if(any(in01)){
    q <- q[in01]
    out[in01] <- pgbetap(q/(1-q), c, d, kappa, tau, scale = 1)
  }
  out
}

#' @rdname GBeta
#' @export
rgbeta <- function(n, c, d, kappa, tau, method = "mixture"){ 
  betap <- rgbetap(n, c, d, kappa, tau, scale = 1, method = method)
  betap / (1 + betap)
}

#' @rdname GBeta
#' @export
qgbeta <- function(p, c, d, kappa, tau){
  stopifnot(all(p >= 0 & p <= 1))
  stopifnot(c > 0, d > 0, tau > 0)
  if(kappa == 0 || tau == 1){
    return(qbeta(p, c, d))
  }
  uq(unuran_gbeta_pinv(c, d, kappa, tau), p)
}
