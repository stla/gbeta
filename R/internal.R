#' @importFrom gsl hyperg_2F1
#' @noRd
Gauss2F1 <- function(a, b, c, x){ 
  if(x >= 0 && x <= 1){ # hyperg_2F1 works fine in this situation
    hyperg_2F1(a, b, c, x)
  }else{ # transform to come down to the first situation
    hyperg_2F1(c - a, b, c, 1 - 1/(1 - x)) / (1 - x)^b
  }
}

#' @importFrom gsl hyperg_2F1 lnpoch
#' @importFrom Runuran unuran.discr.new unuran.new unuran.cont.new ur pinv.new
NULL

pmf_unnormalized <- function(k, c, d, kappa, tau){
  # out <- numeric(length(k))
  # positive <- k >= 0
  # k <- k[positive]
  # out[positive] <- 
  #   if(tau < 1){
  #     exp(k*log(1-tau) - lfactorial(k) + 
  #           lnpoch(c+d-kappa,k) + lnpoch(d,k) - lnpoch(c+d,k)) 
  #   }else{
  #     exp(k*log(1-1/tau) - lfactorial(k) + 
  #           lnpoch(c+d-kappa,k) + lnpoch(c,k) - lnpoch(c+d,k))
  #   }
  # out
  if(tau < 1){
    exp(k*log(1-tau) - lfactorial(k) + 
          lnpoch(c+d-kappa,k) + lnpoch(d,k) - lnpoch(c+d,k)) 
  }else{
    exp(k*log(1-1/tau) - lfactorial(k) + 
          lnpoch(c+d-kappa,k) + lnpoch(c,k) - lnpoch(c+d,k))
  }
}

NormalizingConstant <- function(c, d, kappa, tau){
  if(tau < 1){
    hyperg_2F1(d, c+d-kappa, c+d, 1-tau) 
  }else{
    hyperg_2F1(c, c+d-kappa, c+d, 1-1/tau) 
  }
}

Ksampler <- function(n, c, d, kappa, tau){
  # dist <- unuran.discr.new(
  #   pmf = function(k) pmf_unnormalized(k, c, d, kappa, tau),
  #   lb = 0, ub = Inf, sum = NormalizingConstant(c, d, kappa, tau)
  # )
  # unuran <- unuran.new(dist, method = "dgt") 
  unuran <- unuran.new(unuran.discr.new(
    pmf = function(k) pmf_unnormalized(k, c, d, kappa, tau),
    lb = 0, ub = Inf, sum = NormalizingConstant(c, d, kappa, tau)
  ), method = "dgt")
  ur(unuran, n) 
}

unuran_gbetap_arou <- function(c, d, kappa, tau){
  distr <- unuran.cont.new(
    #cdf = function(q) pgbp(q, c, d, kappa, tau),
    pdf = function(phi) dgbetap(phi, c, d, kappa, tau, log = TRUE),
    dpdf = function(x) (c-1)/x + (kappa-c-d)/(1+x) - kappa*tau/(1+tau*x),
    islog = TRUE,
    lb = 0, ub = Inf
  )
  unuran.new(distr, method = "arou")
}

unuran_gbeta_pinv <- function(c, d, kappa, tau){
  pinv.new(
    cdf = function(x) pgbeta(x, c, d, kappa, tau),
    lb = 0, ub = 1
  )
}
