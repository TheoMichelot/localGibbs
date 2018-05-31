#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _localGibbs_nllkLG_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _localGibbs_rsf(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _localGibbs_scalez(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _localGibbs_simLG_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_localGibbs_nllkLG_rcpp", (DL_FUNC) &_localGibbs_nllkLG_rcpp, 12},
    {"_localGibbs_rsf",         (DL_FUNC) &_localGibbs_rsf,          5},
    {"_localGibbs_scalez",      (DL_FUNC) &_localGibbs_scalez,       5},
    {"_localGibbs_simLG_rcpp",  (DL_FUNC) &_localGibbs_simLG_rcpp,   7},
    {NULL, NULL, 0}
};

void R_init_localGibbs(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}