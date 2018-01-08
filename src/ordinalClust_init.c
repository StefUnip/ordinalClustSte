#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ordinalClust_allej(SEXP, SEXP);
extern SEXP _ordinalClust_compare_vec(SEXP, SEXP);
extern SEXP _ordinalClust_ordiemCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ordinalClust_pej(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ordinalClust_pejp1_ej(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ordinalClust_pejp1_yjej(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ordinalClust_pejp1zj1_ej(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ordinalClust_pejp1zj1_yjej(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ordinalClust_pyj_ej(SEXP, SEXP);
extern SEXP _ordinalClust_unsigned_to_signed(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ordinalClust_allej",              (DL_FUNC) &_ordinalClust_allej,              2},
    {"_ordinalClust_compare_vec",        (DL_FUNC) &_ordinalClust_compare_vec,        2},
    {"_ordinalClust_ordiemCpp",          (DL_FUNC) &_ordinalClust_ordiemCpp,          7},
    {"_ordinalClust_pej",                (DL_FUNC) &_ordinalClust_pej,                6},
    {"_ordinalClust_pejp1_ej",           (DL_FUNC) &_ordinalClust_pejp1_ej,           4},
    {"_ordinalClust_pejp1_yjej",         (DL_FUNC) &_ordinalClust_pejp1_yjej,         5},
    {"_ordinalClust_pejp1zj1_ej",        (DL_FUNC) &_ordinalClust_pejp1zj1_ej,        4},
    {"_ordinalClust_pejp1zj1_yjej",      (DL_FUNC) &_ordinalClust_pejp1zj1_yjej,      5},
    {"_ordinalClust_pyj_ej",             (DL_FUNC) &_ordinalClust_pyj_ej,             2},
    {"_ordinalClust_unsigned_to_signed", (DL_FUNC) &_ordinalClust_unsigned_to_signed, 1},
    {NULL, NULL, 0}
};

void R_init_ordinalClust(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
