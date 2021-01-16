// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rescale_to_coverage
DataFrame rescale_to_coverage(NumericVector cdf, double coverage, NumericVector Mid_Week);
RcppExport SEXP _mypackage_rescale_to_coverage(SEXP cdfSEXP, SEXP coverageSEXP, SEXP Mid_WeekSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type cdf(cdfSEXP);
    Rcpp::traits::input_parameter< double >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Mid_Week(Mid_WeekSEXP);
    rcpp_result_gen = Rcpp::wrap(rescale_to_coverage(cdf, coverage, Mid_Week));
    return rcpp_result_gen;
END_RCPP
}
// incremental_cov_leakage
NumericMatrix incremental_cov_leakage(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _mypackage_incremental_cov_leakage(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(incremental_cov_leakage(x, y));
    return rcpp_result_gen;
END_RCPP
}
// impact_time_week
NumericMatrix impact_time_week(NumericMatrix x, NumericMatrix y, NumericMatrix z, NumericVector a, NumericVector b, NumericVector c);
RcppExport SEXP _mypackage_impact_time_week(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(impact_time_week(x, y, z, a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// cumulative_covC
NumericMatrix cumulative_covC(NumericMatrix x);
RcppExport SEXP _mypackage_cumulative_covC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cumulative_covC(x));
    return rcpp_result_gen;
END_RCPP
}
// incremental_covC
NumericMatrix incremental_covC(NumericMatrix x);
RcppExport SEXP _mypackage_incremental_covC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(incremental_covC(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mypackage_rescale_to_coverage", (DL_FUNC) &_mypackage_rescale_to_coverage, 3},
    {"_mypackage_incremental_cov_leakage", (DL_FUNC) &_mypackage_incremental_cov_leakage, 2},
    {"_mypackage_impact_time_week", (DL_FUNC) &_mypackage_impact_time_week, 6},
    {"_mypackage_cumulative_covC", (DL_FUNC) &_mypackage_cumulative_covC, 1},
    {"_mypackage_incremental_covC", (DL_FUNC) &_mypackage_incremental_covC, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_mypackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
