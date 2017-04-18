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
extern SEXP phylobase_all_naC(SEXP);
extern SEXP phylobase_any_naC(SEXP);
extern SEXP phylobase_checkTreeCpp(SEXP, SEXP);
extern SEXP phylobase_edgeIdCpp(SEXP, SEXP);
extern SEXP phylobase_getAllNodesFast(SEXP);
extern SEXP phylobase_getAllNodesSafe(SEXP);
extern SEXP phylobase_getRange(SEXP, SEXP);
extern SEXP phylobase_hasDuplicatedLabelsCpp(SEXP);
extern SEXP phylobase_hasPolytomy(SEXP);
extern SEXP phylobase_hasSingleton(SEXP);
extern SEXP phylobase_isLabelName(SEXP, SEXP);
extern SEXP phylobase_nb_naC(SEXP);
extern SEXP phylobase_nRoots(SEXP);
extern SEXP phylobase_nTipsFastCpp(SEXP);
extern SEXP phylobase_nTipsSafe(SEXP);
extern SEXP phylobase_tabulateTips(SEXP);
extern SEXP phylobase_testEqInt(SEXP, SEXP);
extern SEXP phylobase_tipsFast(SEXP);
extern SEXP phylobase_tipsSafe(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"phyloxx",       (DL_FUNC) &phyloxx,       6},
    {"reorderBinary", (DL_FUNC) &reorderBinary, 7},
    {"reorderRobust", (DL_FUNC) &reorderRobust, 6},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"ancestors_c",                      (DL_FUNC) &ancestors_c,                      3},
    {"descendants_c",                    (DL_FUNC) &descendants_c,                    3},
    {"phylobase_all_naC",                (DL_FUNC) &phylobase_all_naC,                1},
    {"phylobase_any_naC",                (DL_FUNC) &phylobase_any_naC,                1},
    {"phylobase_checkTreeCpp",           (DL_FUNC) &phylobase_checkTreeCpp,           2},
    {"phylobase_edgeIdCpp",              (DL_FUNC) &phylobase_edgeIdCpp,              2},
    {"phylobase_getAllNodesFast",        (DL_FUNC) &phylobase_getAllNodesFast,        1},
    {"phylobase_getAllNodesSafe",        (DL_FUNC) &phylobase_getAllNodesSafe,        1},
    {"phylobase_getRange",               (DL_FUNC) &phylobase_getRange,               2},
    {"phylobase_hasDuplicatedLabelsCpp", (DL_FUNC) &phylobase_hasDuplicatedLabelsCpp, 1},
    {"phylobase_hasPolytomy",            (DL_FUNC) &phylobase_hasPolytomy,            1},
    {"phylobase_hasSingleton",           (DL_FUNC) &phylobase_hasSingleton,           1},
    {"phylobase_isLabelName",            (DL_FUNC) &phylobase_isLabelName,            2},
    {"phylobase_nb_naC",                 (DL_FUNC) &phylobase_nb_naC,                 1},
    {"phylobase_nRoots",                 (DL_FUNC) &phylobase_nRoots,                 1},
    {"phylobase_nTipsFastCpp",           (DL_FUNC) &phylobase_nTipsFastCpp,           1},
    {"phylobase_nTipsSafe",              (DL_FUNC) &phylobase_nTipsSafe,              1},
    {"phylobase_tabulateTips",           (DL_FUNC) &phylobase_tabulateTips,           1},
    {"phylobase_testEqInt",              (DL_FUNC) &phylobase_testEqInt,              2},
    {"phylobase_tipsFast",               (DL_FUNC) &phylobase_tipsFast,               1},
    {"phylobase_tipsSafe",               (DL_FUNC) &phylobase_tipsSafe,               2},
    {NULL, NULL, 0}
};

void R_init_phylobase(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
