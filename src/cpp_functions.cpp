#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]] 
double det_sympd_C (arma::mat x) {
  arma::mat cholx = chol(x);
  double y = prod(cholx.diag());
  return y * y;
}

// [[Rcpp::export(".dmnorm_C")]]
double dmnorm_C (arma::vec y, arma::vec mu, arma::mat cov_prec, bool is_cov) {
  double n = y.n_elem;
  arma::mat prec;
  double pi = arma::datum::pi;
  if (is_cov){
    prec = inv_sympd(cov_prec);
  } else {
    prec = cov_prec;
  }
  
  arma::mat log_exp = trans(y - mu) * prec * (y - mu);
  double out = -n / 2 * log(2 * pi) + .5 * log(det_sympd_C(prec)) - .5 * log_exp(0, 0);

  return out;
}

// [[Rcpp::export(".update_gp_mean_C")]]
arma::vec update_gp_mean_C(arma::vec y, arma::mat R12, arma::mat R22,
                   arma::vec mu1, arma::vec mu2){
  arma::mat R22i = inv_sympd(R22);
  return mu1 + R12 * (R22i * (y - mu2));
}

// [[Rcpp::export(".update_gp_var_C")]]
arma::mat update_gp_var_C(arma::mat R11, arma::mat R12, arma::mat R22){
  arma::mat R22i = inv_sympd(R22);
  return R11 - R12 * R22i * R12.t();
}
