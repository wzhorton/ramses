#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat cinv(arma::mat x) {
  return inv_sympd(x);
}

