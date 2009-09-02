/*
  descendants.c:
    Identify all descendants of a given node. Function inputs are
  derived from a phylo4 edge matrix, which *must* be in preorder order.
  The isDescendant input vector should contain 1 for the immediate
  children of the node, and 0 otherwise. The function returns this
  vector updated to include all further descendants.
*/

#include <R.h>

void descendants(int *isDescendant, int *ancestor, int *descendant, 
    int *numEdges) {

    int child=0;
    for (int i=0; i<*numEdges; i++) {
        if (isDescendant[i]==1) {
            child = descendant[i];
            for (int j=i+1; j<*numEdges; j++) {
                if (ancestor[j]==child) isDescendant[j]=1; 
            }
        }
    }
}
