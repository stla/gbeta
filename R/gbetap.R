#' @useDynLib gbeta
#' @importFrom Rcpp evalCpp
NULL

#' Generalized Beta prime distribution
#' @description Density, distribution function, quantile function, and random 
#'   generation for the generalized Beta prime distribution.
#'   
#' @param x numeric vector
#' @param q numeric vector of quantiles
#' @param p numeric vector of probabilities
#' @param n positive integer, the desired number of simulations
#' @param c,d,kappa,tau parameters; they must be strictly positive numbers, 
#'   except \code{kappa} which can take any value
#' @param scale scale parameter, a strictly positive number
#' @param log logical, whether to return the log-density
#' @param method the method of random generation, \code{"mixture"} or 
#'   \code{"arou"}; only a positive \code{kappa} is allowed for the 
#'   \code{"mixture"} method, but this method is faster
#'   
#' @references 
#' \itemize{
#'   \item Stéphane Laurent. 
#'     \emph{Some Poisson mixtures distributions with a hyperscale parameter}.
#'     Brazilian Journal of Probability and Statistics 26, No. 3 (2012), 
#'     pp. 265-278.
#'     <doi:10.1214/11-BJPS139>
#'   \item Myriam Chabot.
#'     \emph{Sur l’estimation du rapport de deux paramètres d’intensité
#'           poissonniens et l’estimation par fonctions de masse prédictives}. 
#'     Master thesis. 
#'     Université de Scherbrooke, 2016.
#' }
#'
#' @examples library(gbeta)
#' curve(dgbetap(x, 4, 12, 10, 0.01), to = 10, axes = FALSE, lwd = 2)
#' axis(1)
#'   
#' @rdname GBetaP
#' @name GBetaP
#' @export
dgbetap <- function(x, c, d, kappa, tau, scale = 1, log = FALSE){ # tau rate parameter NO! --> add scale parameter
  stopifnot(c > 0, d > 0, tau > 0)
  out <- numeric(length(x))
  negs <- x < 0
  x <- x[!negs]
  if(length(x)){
    if(kappa == c + d){
      dnsty <- dbetaprime(x, c, d, scale/tau, log)
    }else if(kappa == 0 || tau == 1){
      dnsty <- dbetaprime(x, c, d, scale, log)
    }else{
      logunnormalized <- -lbeta(c,d) + (c-1)*log(x/scale) + 
        (kappa-c-d)*log(1+x/scale) - kappa*log(1+tau*x/scale) - log(scale)
      dnsty <- if(log){
        logunnormalized - log(Gauss2F1(c, kappa, c+d, 1-tau))
      }else{
        exp(logunnormalized) / Gauss2F1(c, kappa, c+d, 1-tau)
      }
      # dnsty <- 1/beta(c,d) * x^(c-1)*(1+x)^(kappa-c-d)/(1+tau*x)^(kappa) / 
      #   Gauss2F1(c, kappa, c+d, 1-tau)
    }
    out[!negs] <- dnsty
  }
  out
}

#' @rdname GBetaP
#' @export
pgbetap <- function(q, c, d, kappa, tau, scale = 1){ # quid if q < 0 ? - done in cpp
  stopifnot(c > 0, d > 0, tau > 0)
  if(kappa == c + d){
    pbetaprime(q, c, d, scale/tau)
  }else if(kappa == 0 || tau == 1){
    pbetaprime(q, c, d, scale)
  }else{
    C <- beta(c,d) * Gauss2F1(c, kappa, c+d, 1-tau)
    (q/scale)^c * euler(c, d, -kappa+c+d, 1/tau, q/scale) / C
  }
}

#' @rdname GBetaP
#' @importFrom Runuran ur
#' @export
rgbetap <- function(n, c, d, kappa, tau, scale = 1, method = "mixture"){
  method <- match.arg(method, c("mixture", "arou"))
  if(method == "mixture" && kappa < 0){
    stop(
      "`kappa < 0` is not allowed with the 'mixture' method."
    )
  }
  stopifnot(c > 0, d > 0, tau > 0)
  if(tau == 1 || kappa == 0){
    rbetaprime(n, c, d, scale)
  }else if(kappa == c + d){
    rbetaprime(n, c, d, scale/tau)
  }else if(method == "mixture"){
    K <- Ksampler(n, c, d, -kappa+c+d, 1/tau)
    if(tau > 1){
      rbetaprime(n, c, d + K, scale)
    }else{
      rbetaprime(n, c + K, d, scale)
    }
  }else if(method == "arou"){
    scale * ur(unuran_gbetap_arou(c, d, kappa, tau), n)
  }
}

#' @rdname GBetaP
#' @export
qgbetap <- function(p, c, d, kappa, tau, scale = 1){
  q <- qgbeta(p, c, d, kappa, tau)
  scale * q / (1 - q)
}
