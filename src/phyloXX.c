/*
  descendants.c:
    Identify all descendants of a given node. Function inputs are
  derived from a phylo4 edge matrix, which *must* be in preorder order.
  The isDescendant input vector should contain 1 for the immediate
  children of the node, and 0 otherwise. The function returns this
  vector updated to include all further descendants.
*/

// test1 <- function() {
//     for (i in edge[, 2]) {
//         dex <- edge[, 1] == i
//         cur <- edge[, 2] == i
//         xx[dex] <- phy@edge.length[dex] + xx[cur]
//         segs$v0x[dex] <- xx[cur]
//     }
//     return(list(segs=segs, xx=xx))
// }
// test1out <- test1()
// segs <- test1out$segs
// xx   <- test1out$xx

// test2 <- function() {
//     for(i in rev((Ntips + 1):nEdges(phy))) {
//         dex <- edge[, 1] == i
//         cur <- edge[, 2] == i
//         yy[cur] <- segs$v0y[dex] <- mean(yy[dex])
//     }
//     return(list(segs=segs, yy=yy))
// }
// test2out <- test2()
// segs <- test2out$segs
// yy   <- test2out$yy
// segs$h0y <- segs$h1y <- segs$v1y <- yy

#include <R.h>

// 
// void phyloyy(int *edge1, int *edge2, int *ntips, 
//                  int *numEdges, double *yy, double *v0y) 
// {
//     int i;
//     int k;
//     int j;
//     int cur;
//     int des;
//     int count;
//     double tmp;
//     double theMean;
//     Rprintf("test\n");
//     for (i=*numEdges; i > *ntips ; i--) {
//         for (k=0; k<*numEdges; k++) {
//             if(i == edge2[k]) {
//                 cur = k;
//             }
//         }
//         tmp=0;
//         count=0;
//         for (j=0; j<*numEdges; j++) {
//             if(i == edge1[j]) {
//                 des = j;
//                 tmp += yy[j];
//                 count += 1;
//             }
//         }
//         theMean  = tmp / count;
//         yy[cur]  = theMean;
//         for (j=0; j<*numEdges; j++) {
//             if(i == edge1[j]) {
//                 v0y[j] = theMean;
//             }
//         }
// 
//     }
// }

void phyloxx(int *edge1, int *edge2, double *edgeLengths, 
                 int *numEdges, double *xx, double *v0x) 
{
    int j;
    int i;
    int k;
    int cur=0;
    for (i=0; i <*numEdges; i++) {
        for (k=0; k<*numEdges; k++) {
            if(edge2[i] == edge2[k]) {
                cur = k;
            }
        }
        for (j=0; j<*numEdges; j++) {
            if(edge2[i] == edge1[j]) {
                xx[j] = edgeLengths[j] + xx[cur];
                v0x[j] = xx[cur];
            }
        }
    }
}
