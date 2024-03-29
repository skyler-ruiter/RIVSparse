// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/RIVSparse.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// convertSparse
void convertSparse(S4& mat);
RcppExport SEXP _RIVSparse_convertSparse(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4& >::type mat(matSEXP);
    convertSparse(mat);
    return R_NilValue;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _RIVSparse_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_vcsc();
RcppExport SEXP _rcpp_module_boot_yada();

static const R_CallMethodDef CallEntries[] = {
    {"_RIVSparse_convertSparse", (DL_FUNC) &_RIVSparse_convertSparse, 1},
    {"_RIVSparse_rcpp_hello_world", (DL_FUNC) &_RIVSparse_rcpp_hello_world, 0},
    {"_rcpp_module_boot_vcsc", (DL_FUNC) &_rcpp_module_boot_vcsc, 0},
    {"_rcpp_module_boot_yada", (DL_FUNC) &_rcpp_module_boot_yada, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_RIVSparse(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
