// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cumz
NumericVector cumz(NumericVector x);
RcppExport SEXP _lemur_cumz(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cumz(x));
    return rcpp_result_gen;
END_RCPP
}
// cumz_which_abs_max
List cumz_which_abs_max(NumericVector x);
RcppExport SEXP _lemur_cumz_which_abs_max(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cumz_which_abs_max(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lemur_cumz", (DL_FUNC) &_lemur_cumz, 1},
    {"_lemur_cumz_which_abs_max", (DL_FUNC) &_lemur_cumz_which_abs_max, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_lemur(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
