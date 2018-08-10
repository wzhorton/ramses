// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// det_sympd_C
double det_sympd_C(arma::mat x);
RcppExport SEXP _ramses_det_sympd_C(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(det_sympd_C(x));
    return rcpp_result_gen;
END_RCPP
}
// dmnorm_C
double dmnorm_C(arma::vec y, arma::vec mu, arma::mat cov_prec, bool is_cov);
RcppExport SEXP _ramses_dmnorm_C(SEXP ySEXP, SEXP muSEXP, SEXP cov_precSEXP, SEXP is_covSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_prec(cov_precSEXP);
    Rcpp::traits::input_parameter< bool >::type is_cov(is_covSEXP);
    rcpp_result_gen = Rcpp::wrap(dmnorm_C(y, mu, cov_prec, is_cov));
    return rcpp_result_gen;
END_RCPP
}
// update_gp_C
arma::vec update_gp_C(arma::vec y, arma::mat R12, arma::mat R22, arma::vec mu1, arma::vec mu2);
RcppExport SEXP _ramses_update_gp_C(SEXP ySEXP, SEXP R12SEXP, SEXP R22SEXP, SEXP mu1SEXP, SEXP mu2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R12(R12SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R22(R22SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu2(mu2SEXP);
    rcpp_result_gen = Rcpp::wrap(update_gp_C(y, R12, R22, mu1, mu2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ramses_det_sympd_C", (DL_FUNC) &_ramses_det_sympd_C, 1},
    {"_ramses_dmnorm_C", (DL_FUNC) &_ramses_dmnorm_C, 4},
    {"_ramses_update_gp_C", (DL_FUNC) &_ramses_update_gp_C, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_ramses(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
