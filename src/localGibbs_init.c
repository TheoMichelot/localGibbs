#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Generated with tools::package_native_routine_registration_skeleton

/* .Call calls */
extern SEXP _localGibbs_nllkLG_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _localGibbs_rsf(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _localGibbs_scalez(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _localGibbs_simLG_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _localGibbs_simSSF_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _localGibbs_simZeros_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_localGibbs_nllkLG_rcpp",   (DL_FUNC) &_localGibbs_nllkLG_rcpp,   12},
    {"_localGibbs_rsf",           (DL_FUNC) &_localGibbs_rsf,            5},
    {"_localGibbs_scalez",        (DL_FUNC) &_localGibbs_scalez,         5},
    {"_localGibbs_simLG_rcpp",    (DL_FUNC) &_localGibbs_simLG_rcpp,     7},
    {"_localGibbs_simSSF_rcpp",   (DL_FUNC) &_localGibbs_simSSF_rcpp,   12},
    {"_localGibbs_simZeros_rcpp", (DL_FUNC) &_localGibbs_simZeros_rcpp,  6},
    {NULL, NULL, 0}
};

void R_init_localGibbs(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}