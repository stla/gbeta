// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;
#include <Rcpp.h>

//---------- integrand ---------//
double integrand(double u,
                 double c,
                 double d,
                 double kappa,
                 double tau,
                 double q) {
  double f = pow(u, c - 1) * pow(1 + q * u, -kappa) *
             pow(1 + q * u / tau, kappa - c - d);
  return f;
}

class Integrand : public Func {
 private:
  double c;
  double d;
  double kappa;
  double tau;
  double q;

 public:
  Integrand(double c_, double d_, double kappa_, double tau_, double q_)
      : c(c_), d(d_), kappa(kappa_), tau(tau_), q(q_) {}

  double operator()(const double& u) const {
    return integrand(u, c, d, kappa, tau, q);
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector euler(const double c,
                          const double d,
                          const double kappa,
                          const double tau,
                          const Rcpp::NumericVector q,
                          const int subdiv = 100,
                          const double eps_abs = 1e-14,
                          const double eps_rel = 1e-14) {
  const size_t n = q.size();
  Rcpp::NumericVector err_ests(n);
  Rcpp::IntegerVector err_codes(n);
  Rcpp::NumericVector out = Rcpp::NumericVector(n);
  for(size_t i = 0; i < n; i++) {
    if(q(i) <= 0){
      out(i) = 0;
      err_ests(i) = 0;
      err_codes(i) = 0;
    }else{
      const Integrand f(c, d, kappa, tau, q(i));
      double err_est;
      int err_code;
      const double res = integrate(f, 0, 1, err_est, err_code, subdiv, eps_abs,
                                   eps_rel, Integrator<double>::GaussKronrod201);
      out(i) = res;
      err_ests(i) = err_est;
      err_codes(i) = err_code;
    }
  }
  out.attr("err_est") = err_ests;
  out.attr("err_code") = err_codes;
  return out;
}