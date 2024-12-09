// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// hmm_em
List hmm_em(NumericMatrix observations, int n_states, std::string dist_type, int max_iter);
RcppExport SEXP _wpfSA24204165_hmm_em(SEXP observationsSEXP, SEXP n_statesSEXP, SEXP dist_typeSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< int >::type n_states(n_statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_type(dist_typeSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(hmm_em(observations, n_states, dist_type, max_iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_wpfSA24204165_hmm_em", (DL_FUNC) &_wpfSA24204165_hmm_em, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_wpfSA24204165(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}