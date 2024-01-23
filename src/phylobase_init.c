#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void phyloxx(void *, void *, void *, void *, void *, void *);
extern void reorderBinary(void *, void *, void *, void *, void *, void *, void *);
extern void reorderRobust(void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP ancestors_c(SEXP, SEXP, SEXP);
extern SEXP descendants_c(SEXP, SEXP, SEXP);
extern SEXP _phylobase_all_naC(SEXP);
extern SEXP _phylobase_any_naC(SEXP);
extern SEXP _phylobase_checkTreeCpp(SEXP, SEXP);
extern SEXP _phylobase_edgeIdCpp(SEXP, SEXP);
extern SEXP _phylobase_getAllNodesFast(SEXP);
extern SEXP _phylobase_getAllNodesSafe(SEXP);
extern SEXP _phylobase_getRange(SEXP, SEXP);
extern SEXP _phylobase_hasDuplicatedLabelsCpp(SEXP);
extern SEXP _phylobase_hasPolytomy(SEXP);
extern SEXP _phylobase_hasSingleton(SEXP);
extern SEXP _phylobase_isLabelName(SEXP, SEXP);
extern SEXP _phylobase_nb_naC(SEXP);
extern SEXP _phylobase_nRoots(SEXP);
extern SEXP _phylobase_nTipsFastCpp(SEXP);
extern SEXP _phylobase_nTipsSafe(SEXP);
extern SEXP _phylobase_tabulateTips(SEXP);
extern SEXP _phylobase_testEqInt(SEXP, SEXP);
extern SEXP _phylobase_tipsFast(SEXP);
extern SEXP _phylobase_tipsSafe(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"phyloxx",       (DL_FUNC) &phyloxx,       6},
    {"reorderBinary", (DL_FUNC) &reorderBinary, 7},
    {"reorderRobust", (DL_FUNC) &reorderRobust, 6},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"ancestors_c",                      (DL_FUNC) &ancestors_c,                      3},
    {"descendants_c",                    (DL_FUNC) &descendants_c,                    3},
    {"_phylobase_all_naC",                (DL_FUNC) &_phylobase_all_naC,                1},
    {"_phylobase_any_naC",                (DL_FUNC) &_phylobase_any_naC,                1},
    {"_phylobase_checkTreeCpp",           (DL_FUNC) &_phylobase_checkTreeCpp,           2},
    {"_phylobase_edgeIdCpp",              (DL_FUNC) &_phylobase_edgeIdCpp,              2},
    {"_phylobase_getAllNodesFast",        (DL_FUNC) &_phylobase_getAllNodesFast,        1},
    {"_phylobase_getAllNodesSafe",        (DL_FUNC) &_phylobase_getAllNodesSafe,        1},
    {"_phylobase_getRange",               (DL_FUNC) &_phylobase_getRange,               2},
    {"_phylobase_hasDuplicatedLabelsCpp", (DL_FUNC) &_phylobase_hasDuplicatedLabelsCpp, 1},
    {"_phylobase_hasPolytomy",            (DL_FUNC) &_phylobase_hasPolytomy,            1},
    {"_phylobase_hasSingleton",           (DL_FUNC) &_phylobase_hasSingleton,           1},
    {"_phylobase_isLabelName",            (DL_FUNC) &_phylobase_isLabelName,            2},
    {"_phylobase_nb_naC",                 (DL_FUNC) &_phylobase_nb_naC,                 1},
    {"_phylobase_nRoots",                 (DL_FUNC) &_phylobase_nRoots,                 1},
    {"_phylobase_nTipsFastCpp",           (DL_FUNC) &_phylobase_nTipsFastCpp,           1},
    {"_phylobase_nTipsSafe",              (DL_FUNC) &_phylobase_nTipsSafe,              1},
    {"_phylobase_tabulateTips",           (DL_FUNC) &_phylobase_tabulateTips,           1},
    {"_phylobase_testEqInt",              (DL_FUNC) &_phylobase_testEqInt,              2},
    {"_phylobase_tipsFast",               (DL_FUNC) &_phylobase_tipsFast,               1},
    {"_phylobase_tipsSafe",               (DL_FUNC) &_phylobase_tipsSafe,               2},
    {NULL, NULL, 0}
};

void R_init_phylobase(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
