// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <Rcpp.h>
#include <algorithm>    // std::count_if
#include <vector>       // std::vector
#include <string>       //

template <typename T>
std::string NumberToString ( T Number ) {
    std::ostringstream ss;
    ss << Number;
    return ss.str();
}

bool isZero(int i) { return (i == 0); }
bool isOne(int i) { return ( i == 1); }
bool isSupTwo(int i) { return (i > 2); }
bool isEqual(int i, int j) { return (i == j); }

Rcpp::IntegerVector getAnces(Rcpp::IntegerMatrix obj) {
// returns the first column (ancestors) of the edge matrix
    Rcpp::IntegerMatrix::Column out = obj( Rcpp::_ , 0);
    return out;
}

Rcpp::IntegerVector getDesc(Rcpp::IntegerMatrix obj) {
// returns the second column (descendants) of the edge matrix
    Rcpp::IntegerMatrix::Column out = obj( Rcpp::_ , 1);
    return out;
}

//[[Rcpp::export]]
bool isLabelName(Rcpp::CharacterVector lblToCheck,
		 Rcpp::CharacterVector lbl ) {

    Rcpp::CharacterVector noLbl = Rcpp::setdiff(lblToCheck, lbl);
    return noLbl.size() == 0;
}

//[[Rcpp::export]]
int nRoots (Rcpp::IntegerVector ances) {
    int ans = std::count (ances.begin(), ances.end(), 0);
    return ans;
}

//[[Rcpp::export]]
std::vector<int> tabulateTips (Rcpp::IntegerVector ances) {
// tabulates ancestor nodes that are not the root.
    int n = Rcpp::max(ances);
    std::vector<int> ans(n);
    for (int i=0; i < ances.size(); i++) {
        int j = ances[i];
        if (j > 0) {
            ans[j - 1]++;
        }
    }
    return ans;
}

//[[Rcpp::export]]
int nTipsSafe (Rcpp::IntegerVector ances) {
// count how many zeros are in the tabulated vector of ancestors
// this gives the number of tips
    std::vector<int> tabTips = tabulateTips(ances);
    int j = count_if (tabTips.begin(), tabTips.end(), isZero);
    return j;
}

//[[Rcpp::export]]
int nTipsFastCpp (Rcpp::IntegerVector ances) {
// if nodes are correctly numbered min(ances) - 1 = nb of tips
// (after removing the root, which is equal to 0).
    int nroots = nRoots(ances);
    if (nroots > 0) {
	int whichRoot = Rcpp::which_min(ances);
	ances.erase(whichRoot);
    }
    int tmp = Rcpp::min(ances);
    return tmp - 1;
}

//[[Rcpp::export]]
bool hasSingleton (Rcpp::IntegerVector ances) {
    std::vector<int> tabTips = tabulateTips(ances);
    int j = count_if (tabTips.begin(), tabTips.end(), isOne);
    return j > 0;
}

//[[Rcpp::export]]
bool hasPolytomy (Rcpp::IntegerVector ances) {
    std::vector<int> tabTips = tabulateTips(ances);
    int j = count_if (tabTips.begin(), tabTips.end(), isSupTwo);
    return j > 0;
}


//[[Rcpp::export]]
Rcpp::IntegerVector tipsSafe (Rcpp::IntegerVector ances, Rcpp::IntegerVector desc) {
    Rcpp::IntegerVector res = Rcpp::match(desc, ances);
    Rcpp::LogicalVector istip = Rcpp::is_na(res);
    int nedge = ances.size();
    std::vector<int> y(nedge);
    int j = 0;
    for(int i = 0; i < nedge; i++) {
	if (istip[i]) {
	    y[j] = desc[i];
	    j++;
	}
    }
    Rcpp::IntegerVector ans(j);
    std::copy (y.begin(), y.begin()+j, ans.begin());
    std::sort  (ans.begin(), ans.end());
    return ans;
}

//[[Rcpp::export]]
Rcpp::IntegerVector tipsFast (Rcpp::IntegerVector ances) {
    int ntips = nTipsFastCpp(ances);
    Rcpp::IntegerVector ans = Rcpp::seq_len(ntips);
    return ans;
}


//[[Rcpp::export]]
Rcpp::IntegerVector getAllNodesSafe (Rcpp::IntegerMatrix edge) {
    Rcpp::IntegerVector ans = Rcpp::as_vector(edge);
    Rcpp::IntegerVector tmp = Rcpp::unique(ans);
    std::sort(tmp.begin(), tmp.end());
    return tmp;
}

//[[Rcpp::export]]
Rcpp::IntegerVector getAllNodesFast (Rcpp::IntegerMatrix edge) {
    Rcpp::IntegerVector tmp = Rcpp::as_vector(edge);
    Rcpp::IntegerVector maxN = Rcpp::range(tmp);
    Rcpp::IntegerVector ans;
    if (maxN[0] == 0) {
        ans = Rcpp::seq_len(maxN[1] + 1);
        ans = ans - 1;
    }
    else {
        ans = Rcpp::seq_len(maxN[1]);
    }
    return ans;
}


// Rcpp::List testNodes (Rcpp::IntegerMatrix edge, bool rooted) {
//     Rcpp::IntegerVector allNodes = Rcpp::as_vector(edge);
//     allNodes = Rcpp::unique(allNodes);
//     std::sort (allNodes.begin(), allNodes.end());
//     Rcpp::IntegerVector supposedNodes = getAllNodesFast(edge, rooted);
//     Rcpp::IntegerVector test = Rcpp::setdiff(supposedNodes, allNodes);
//     Rcpp::LogicalVector res = supposedNodes == allNodes;
//     return Rcpp::List::create(supposedNodes, allNodes, test, res);
// }

//[[Rcpp::export]]
Rcpp::List testEqInt (Rcpp::IntegerVector x, Rcpp::IntegerVector y) {
    Rcpp::LogicalVector xy = x == y;
    Rcpp::LogicalVector yx = y == x;
    return Rcpp::List::create(xy, yx);
}

// Rcpp::IntegerVector getInternalNodes (Rcpp::IntegerMatrix edge, bool rooted) {
//     Rcpp::IntegerVector ances = getAnces(edge);
//     Rcpp::IntegerVector allNodes = getAllNodesFast(edge, rooted);
//     Rcpp::IntegerVector tips = tipsFast(ances);
//     Rcpp::IntegerVector intNodes = Rcpp::setdiff(allNodes, tips);
//     intNodes.erase(intNodes.begin());
//     return intNodes;
// }

//[[Rcpp::export]]
bool all_naC (Rcpp::NumericVector x) {
    return is_true(all(is_na(x)));
}

//[[Rcpp::export]]
bool any_naC (Rcpp::NumericVector x) {
    return is_true(any(is_na(x)));
}

//[[Rcpp::export]]
int nb_naC (Rcpp::NumericVector x) {
    return sum(is_na(x));
}


//[[Rcpp::export]]
Rcpp::NumericVector getRange(Rcpp::NumericVector x, const bool na_rm) {
    Rcpp::NumericVector out(2);
    out[0] = R_PosInf;
    out[1] = R_NegInf;

    int n = x.length();
    for(int i = 0; i < n; ++i) {
	if (!na_rm && R_IsNA(x[i])) {
	    out[0] = NA_REAL;
	    out[1] = NA_REAL;
	    return(out);
	}

	if (x[i] < out[0]) out[0] = x[i];
	if (x[i] > out[1]) out[1] = x[i];
    }

    return(out);
}

//[[Rcpp::export]]
bool hasDuplicatedLabelsCpp (Rcpp::CharacterVector label) {
    return is_true(any(Rcpp::duplicated(na_omit(label))));
}

Rcpp::CharacterVector edgeIdCppInternal (Rcpp::IntegerVector tmp1, Rcpp::IntegerVector tmp2) {
    Rcpp::CharacterVector tmpV1 = Rcpp::as< Rcpp::CharacterVector >(tmp1);
    Rcpp::CharacterVector tmpV2 = Rcpp::as< Rcpp::CharacterVector >(tmp2);
    int Ne = tmp1.size();
    Rcpp::CharacterVector res(Ne);
    for (int i = 0; i < Ne; i++) {
        std::string tmpS1;
        tmpS1 = tmpV1[i];
        std::string tmpS2;
        tmpS2 = tmpV2[i];
        std::string tmpS;
        tmpS = tmpS1.append("-");
        tmpS = tmpS.append(tmpS2);
        res[i] = tmpS;
    }
    return res;
}

//[[Rcpp::export]]
Rcpp::CharacterVector edgeIdCpp (Rcpp::IntegerMatrix edge, std::string type) {
    Rcpp::IntegerVector ances = getAnces(edge);
    Rcpp::IntegerVector desc = getDesc(edge);
    int nedge;

    if (type == "tip" || type == "internal") {
	Rcpp::IntegerVector tips = tipsFast(ances);
	nedge = tips.size();
	Rcpp::IntegerVector ans = match(tips, desc);
	if (type == "tip") {
            Rcpp::IntegerVector tmpAnces(nedge);
            Rcpp::IntegerVector tmpDesc(nedge);
	    for (int j = 0; j < nedge; j++) {
                tmpAnces[j] = ances[ans[j]-1];
                tmpDesc[j] = desc[ans[j]-1];
	    }
            Rcpp::CharacterVector c1(nedge);
            c1 = edgeIdCppInternal(tmpAnces, tmpDesc);
            return c1;
	}
	else if (type == "internal") {
	    int allEdges = ances.size();
	    Rcpp::IntegerVector idEdge = Rcpp::seq_len(allEdges);
	    Rcpp::IntegerVector intnd = Rcpp::setdiff(idEdge, ans);
	    nedge = intnd.size();
            Rcpp::IntegerVector tmpAnces(nedge);
            Rcpp::IntegerVector tmpDesc(nedge);
	    for (int j = 0; j < nedge; j++) {
                tmpAnces[j] = ances[intnd[j]-1];
                tmpDesc[j] = desc[intnd[j]-1];
            }
            Rcpp::CharacterVector c1(nedge);
            c1 = edgeIdCppInternal(tmpAnces, tmpDesc);
            return c1;
	}
    }
    else {
        nedge = ances.size();
        Rcpp::IntegerVector tmpAnces = ances;
        Rcpp::IntegerVector tmpDesc = desc;
        Rcpp::CharacterVector c1(nedge);
        c1 = edgeIdCppInternal(tmpAnces, tmpDesc);
        return c1;
    }
    return "";
}

//[[Rcpp::export]]
Rcpp::List checkTreeCpp(Rcpp::S4 obj, Rcpp::List opts) {

    std::string err, wrn;
    Rcpp::IntegerMatrix ed = obj.slot("edge");
    int nrow = ed.nrow();
    Rcpp::IntegerVector ances = getAnces(ed);
    //Rcpp::IntegerVector desc = getDesc(ed);
    int nroots = nRoots(ances);
    //bool rooted = nroots > 0;
    Rcpp::NumericVector edLength = obj.slot("edge.length");
    Rcpp::CharacterVector edLengthNm = edLength.names();
    Rcpp::CharacterVector label = obj.slot("label");
    Rcpp::CharacterVector labelNm = label.names();
    Rcpp::CharacterVector edLabel = obj.slot("edge.label");
    Rcpp::CharacterVector edLabelNm = edLabel.names();
    Rcpp::IntegerVector allnodesSafe = getAllNodesSafe(ed);
    Rcpp::IntegerVector allnodesFast = getAllNodesFast(ed);
    int nEdLength = edLength.size();
    //int nLabel = label.size();
    //int nEdLabel = edLabel.size();
    int nEdges = nrow;
    bool hasEdgeLength = !all_naC(edLength);

    // check tips
    int ntipsSafe = nTipsSafe(ances);
    int ntipsFast = nTipsFastCpp(ances);
    bool testnTips = ntipsFast == ntipsSafe;
    if (! testnTips) {
	err.append("Tips incorrectly labeled. ");
    }

    //check internal nodes
    bool testNodes = Rcpp::all(allnodesSafe == allnodesFast).is_true() && // is both ways comparison needed?
    	Rcpp::all(allnodesFast == allnodesSafe).is_true();
    if (! testNodes) {
    	err.append("Nodes incorrectly labeled. ");
    }

    // check edge lengths
    if (hasEdgeLength) {
    	if (nEdLength != nEdges) {
    	    err.append("Number of edge lengths do not match number of edges. ");
    	}
    	// if (nb_naC(edLength) > nroots) { // not enough!  -- best done in R
    	//     err.append("Only the root should have NA as an edge length. ");
    	// }
    	if (getRange(edLength, TRUE)[0] < 0) {
    	    err.append("Edge lengths must be non-negative. ");
    	}
    	Rcpp::CharacterVector edgeLblSupp = edgeIdCpp(ed, "all");
	Rcpp::CharacterVector edgeLblDiff = Rcpp::setdiff(edLengthNm, edgeLblSupp);
    	if ( edgeLblDiff.size() != 0 ) {
    	    err.append("Edge lengths incorrectly labeled. ");
    	}
    }

    // check label names
    Rcpp::CharacterVector chrLabelNm = Rcpp::as<Rcpp::CharacterVector>(allnodesFast);
    int j = 0;
    while (j < nroots) { //remove root(s)
    	chrLabelNm.erase(0);
    	j++;
    }
    bool testLabelNm = isLabelName(labelNm, chrLabelNm);
    if (!testLabelNm) {
    	err.append("Tip and node labels must be a named vector, the names must match the node IDs. ");
    	err.append("Use tipLabels<- and/or nodeLabels<- to update them. ");
    }

    // check that tips have labels
    Rcpp::CharacterVector tiplabel(ntipsFast);
    std::copy (label.begin(), label.begin()+ntipsFast, tiplabel.begin());
    bool emptyTipLabel = is_true(any(Rcpp::is_na(tiplabel)));
    if ( emptyTipLabel ) {
    	err.append("All tips must have a label.");
    }

    // check edgeLabels
    Rcpp::CharacterVector chrEdgeLblNm = edgeIdCpp(ed, "all");
    bool testEdgeLblNm = isLabelName(edLabelNm, chrEdgeLblNm);
    if (!testEdgeLblNm) {
    	err.append("Edge labels are not labelled correctly. Use the function edgeLabels<- to update them. ");
    }

    // make sure that tips and node labels are unique
    if (hasDuplicatedLabelsCpp(label)) {
	std::string labOpt = opts["allow.duplicated.labels"];
	if (labOpt == "fail") {
	    err.append("Labels are not unique. ");
	}
	if (labOpt == "warn") {
	    wrn.append("Labels are not unique. ");
	}
    }

    // check for polytomies
    if (hasPolytomy(ances)) {
	std::string msgPoly = "Tree includes polytomies. ";
	std::string polyOpt = opts["poly"];
	if (polyOpt == "fail") {
	    err.append(msgPoly);
	}
	if (polyOpt == "warn") {
	    wrn.append(msgPoly);
	}
    }

    // check number of roots
    if (nroots > 1) {
	std::string msgRoot = "Tree has more than one root. ";
	std::string rootOpt = opts["multiroot"];
	if (rootOpt == "fail") {
	    err.append(msgRoot);
	}
	if (rootOpt == "warn") {
	    wrn.append(msgRoot);
	}
    }

    // check for singletons
    if (hasSingleton(ances)) {
	std::string msgSing = "Tree contains singleton nodes. ";
	std::string singOpt = opts["singleton"];
	if (singOpt == "fail") {
	    err.append(msgSing);
	}
	if (singOpt == "warn") {
	    wrn.append(msgSing);
	}
    }

    return Rcpp::List::create(err, wrn);
}
