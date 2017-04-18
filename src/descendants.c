/*
  descendants.c:
    Identify all descendants of each node in the input vector. Function
  inputs are derived from a phylo4 edge matrix, which *must* be in
  preorder order. The isDescendant output is an indicator matrix of
  which nodes (rows, corresponding to the decendant vector) are
  descendants of each input node (columns, corresponding to the nodes
  vector). It will contain 1 for each descendant of the node, *including
  itself*, and 0 for all other nodes.

  Jim Regetz (NCEAS)
*/

#include <R.h>
#include <Rinternals.h>

SEXP descendants_c(SEXP nod, SEXP anc, SEXP des) {

    int numEdges = length(anc);
    int numNodes = length(nod);

    int* nodes = INTEGER(nod);
    int* ancestor = INTEGER(anc);
    int* descendant = INTEGER(des);

    int child = 0;
    SEXP isDescendant;

    PROTECT(isDescendant = allocMatrix(INTSXP, numEdges, numNodes));
    for (int n=0; n<numNodes; n++) {
        for (int i=0; i<numEdges; i++) {
            if (nodes[n]==descendant[i]) {
                INTEGER(isDescendant)[i + n*numEdges] = 1;
            } else {
                INTEGER(isDescendant)[i + n*numEdges] = 0;
            }
        }
    }
    for (int n=0; n<numNodes; n++) {
        for (int i=0; i<numEdges; i++) {
            if (INTEGER(isDescendant)[i + n*numEdges]==1) {
                child = descendant[i];
                for (int j=i+1; j<numEdges; j++) {
                    if (ancestor[j]==child) {
                        INTEGER(isDescendant)[j + n*numEdges] = 1;
                    }
                }
            }
        }
    }
    UNPROTECT(1);
    return isDescendant;
}
