// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_boottei_parallel
Rcpp::NumericMatrix rcpp_boottei_parallel(const arma::sp_mat& expression, Rcpp::NumericVector ps, const int& permutations, int ncores);
RcppExport SEXP _scTEI_rcpp_boottei_parallel(SEXP expressionSEXP, SEXP psSEXP, SEXP permutationsSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type expression(expressionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< const int& >::type permutations(permutationsSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_boottei_parallel(expression, ps, permutations, ncores));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_meanMatrix_parallel
Rcpp::NumericMatrix rcpp_meanMatrix_parallel(const arma::sp_mat& expression, Rcpp::NumericVector ps, Rcpp::NumericVector psgroup, int ncores);
RcppExport SEXP _scTEI_rcpp_meanMatrix_parallel(SEXP expressionSEXP, SEXP psSEXP, SEXP psgroupSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type expression(expressionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type psgroup(psgroupSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_meanMatrix_parallel(expression, ps, psgroup, ncores));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_pMatrix_parallel
Rcpp::NumericMatrix rcpp_pMatrix_parallel(const arma::sp_mat& expression, Rcpp::NumericVector ps, int ncores);
RcppExport SEXP _scTEI_rcpp_pMatrix_parallel(SEXP expressionSEXP, SEXP psSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type expression(expressionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_pMatrix_parallel(expression, ps, ncores));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_pStrata_parallel
Rcpp::NumericMatrix rcpp_pStrata_parallel(const arma::sp_mat& expression, Rcpp::NumericVector ps, Rcpp::NumericVector psgroup, int ncores);
RcppExport SEXP _scTEI_rcpp_pStrata_parallel(SEXP expressionSEXP, SEXP psSEXP, SEXP psgroupSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type expression(expressionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type psgroup(psgroupSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_pStrata_parallel(expression, ps, psgroup, ncores));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_tei_parallel
Rcpp::List rcpp_tei_parallel(const arma::sp_mat& expression, Rcpp::NumericVector ps, int ncores);
RcppExport SEXP _scTEI_rcpp_tei_parallel(SEXP expressionSEXP, SEXP psSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type expression(expressionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_tei_parallel(expression, ps, ncores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scTEI_rcpp_boottei_parallel", (DL_FUNC) &_scTEI_rcpp_boottei_parallel, 4},
    {"_scTEI_rcpp_meanMatrix_parallel", (DL_FUNC) &_scTEI_rcpp_meanMatrix_parallel, 4},
    {"_scTEI_rcpp_pMatrix_parallel", (DL_FUNC) &_scTEI_rcpp_pMatrix_parallel, 3},
    {"_scTEI_rcpp_pStrata_parallel", (DL_FUNC) &_scTEI_rcpp_pStrata_parallel, 4},
    {"_scTEI_rcpp_tei_parallel", (DL_FUNC) &_scTEI_rcpp_tei_parallel, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_scTEI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
