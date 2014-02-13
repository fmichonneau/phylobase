// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// isLabelName
bool isLabelName(Rcpp::CharacterVector lblToCheck, Rcpp::CharacterVector lbl);
RcppExport SEXP phylobase_isLabelName(SEXP lblToCheckSEXP, SEXP lblSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type lblToCheck(lblToCheckSEXP );
        Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type lbl(lblSEXP );
        bool __result = isLabelName(lblToCheck, lbl);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// nRoots
int nRoots(Rcpp::IntegerVector ances);
RcppExport SEXP phylobase_nRoots(SEXP ancesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP );
        int __result = nRoots(ances);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// tabulateTips
std::vector<int> tabulateTips(Rcpp::IntegerVector ances);
RcppExport SEXP phylobase_tabulateTips(SEXP ancesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP );
        std::vector<int> __result = tabulateTips(ances);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// nTipsSafe
int nTipsSafe(Rcpp::IntegerVector ances);
RcppExport SEXP phylobase_nTipsSafe(SEXP ancesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP );
        int __result = nTipsSafe(ances);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// nTipsFastCpp
int nTipsFastCpp(Rcpp::IntegerVector ances);
RcppExport SEXP phylobase_nTipsFastCpp(SEXP ancesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP );
        int __result = nTipsFastCpp(ances);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// hasSingleton
bool hasSingleton(Rcpp::IntegerVector ances);
RcppExport SEXP phylobase_hasSingleton(SEXP ancesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP );
        bool __result = hasSingleton(ances);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// hasPolytomy
bool hasPolytomy(Rcpp::IntegerVector ances);
RcppExport SEXP phylobase_hasPolytomy(SEXP ancesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP );
        bool __result = hasPolytomy(ances);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// tipsSafe
Rcpp::IntegerVector tipsSafe(Rcpp::IntegerVector ances, Rcpp::IntegerVector desc);
RcppExport SEXP phylobase_tipsSafe(SEXP ancesSEXP, SEXP descSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP );
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type desc(descSEXP );
        Rcpp::IntegerVector __result = tipsSafe(ances, desc);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// tipsFast
Rcpp::IntegerVector tipsFast(Rcpp::IntegerVector ances);
RcppExport SEXP phylobase_tipsFast(SEXP ancesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP );
        Rcpp::IntegerVector __result = tipsFast(ances);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// getAllNodesSafe
Rcpp::IntegerVector getAllNodesSafe(Rcpp::IntegerMatrix edge);
RcppExport SEXP phylobase_getAllNodesSafe(SEXP edgeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type edge(edgeSEXP );
        Rcpp::IntegerVector __result = getAllNodesSafe(edge);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// getAllNodesFast
Rcpp::IntegerVector getAllNodesFast(Rcpp::IntegerMatrix edge, bool rooted);
RcppExport SEXP phylobase_getAllNodesFast(SEXP edgeSEXP, SEXP rootedSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type edge(edgeSEXP );
        Rcpp::traits::input_parameter< bool >::type rooted(rootedSEXP );
        Rcpp::IntegerVector __result = getAllNodesFast(edge, rooted);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// testEqInt
Rcpp::List testEqInt(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
RcppExport SEXP phylobase_testEqInt(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type y(ySEXP );
        Rcpp::List __result = testEqInt(x, y);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// all_naC
bool all_naC(Rcpp::NumericVector x);
RcppExport SEXP phylobase_all_naC(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP );
        bool __result = all_naC(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// any_naC
bool any_naC(Rcpp::NumericVector x);
RcppExport SEXP phylobase_any_naC(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP );
        bool __result = any_naC(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// nb_naC
int nb_naC(Rcpp::NumericVector x);
RcppExport SEXP phylobase_nb_naC(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP );
        int __result = nb_naC(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// getRange
Rcpp::NumericVector getRange(Rcpp::NumericVector x, const bool na_rm);
RcppExport SEXP phylobase_getRange(SEXP xSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< const bool >::type na_rm(na_rmSEXP );
        Rcpp::NumericVector __result = getRange(x, na_rm);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// hasDuplicatedLabelsCpp
bool hasDuplicatedLabelsCpp(Rcpp::CharacterVector label);
RcppExport SEXP phylobase_hasDuplicatedLabelsCpp(SEXP labelSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type label(labelSEXP );
        bool __result = hasDuplicatedLabelsCpp(label);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// edgeIdCpp
Rcpp::CharacterVector edgeIdCpp(Rcpp::IntegerMatrix edge, std::string type);
RcppExport SEXP phylobase_edgeIdCpp(SEXP edgeSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type edge(edgeSEXP );
        Rcpp::traits::input_parameter< std::string >::type type(typeSEXP );
        Rcpp::CharacterVector __result = edgeIdCpp(edge, type);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// checkTreeCpp
Rcpp::List checkTreeCpp(Rcpp::S4 obj, Rcpp::List opts);
RcppExport SEXP phylobase_checkTreeCpp(SEXP objSEXP, SEXP optsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::S4 >::type obj(objSEXP );
        Rcpp::traits::input_parameter< Rcpp::List >::type opts(optsSEXP );
        Rcpp::List __result = checkTreeCpp(obj, opts);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
