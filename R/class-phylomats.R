
#' matrix classes for phylobase
#' 
#' Classes representing phylogenies as matrices
#' 
#' 
#' @name phylomat-class
#' @aliases phylo4vcov-class as_phylo4vcov
#' @docType class
#' @param from a \code{phylo4} object
#' @param \dots optional arguments, to be passed to \code{vcov.phylo} in
#' \code{ape} (the main useful option is \code{cor}, which can be set to
#' \code{TRUE} to compute a correlation rather than a variance-covariance
#' matrix)
#' @section Objects from the Class: These are square matrices (with rows and
#' columns corresponding to tips, and internal nodes implicit) with different
#' meanings depending on the type (variance-covariance matrix, distance matrix,
#' etc.).
#' @author Ben Bolker
#' @keywords classes
#' @examples
#' 
#'   tree.owls <- ape::read.tree(text="(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);")
#'   o2 <- as(tree.owls,"phylo4")
#'   ov <- as(o2,"phylo4vcov")
#'   o3 <- as(ov,"phylo4")
#'   ## these are not completely identical, but are
#'   ## topologically identical ...
#' 
#'   ## edge matrices are in a different order:
#'   ## cf. o2@edge and o3@edge
#'   ## BUT the edge matrices are otherwise identical
#'   identical(o2@edge[order(o2@edge[,2]),],
#'             o3@edge[order(o3@edge[,2]),])
#' 
#'   ## There is left/right ambiguity here in the tree orders:
#'   ## in o2 the 5->6->7->1 lineage
#'   ## (terminating in Strix aluco)
#'   ## is first, in o3 the 5->6->3 lineage
#'   ## (terminating in Athene noctua) is first.
#' 
#' 
## define class for phylogenetic var-cov matrices
setClass("phylo4vcov",
         representation("matrix",
                        edge.label="character",
                        order="character"))

## phylo4 -> var-cov: simply wrap ape::vcv.phylo
##  and add other slots
as_phylo4vcov <- function(from,...) {
  m <- ape::vcv.phylo(as(from,"phylo"),...)
  new("phylo4vcov",
      m,
      edge.label=from@edge.label,
      order=from@order)
}
setAs("phylo4","phylo4vcov",
      function(from,to) {
        as_phylo4vcov(from)})

## var-cov to phylo4
setAs("phylo4vcov","phylo4",
      function(from,to) {
        matrix2tree <- function(v,reorder=TRUE) {
          ## no polytomies allowed
          va <- v
          tipnames <- rownames(v)
          ntip <- nrow(v)
          dimnames(v) <- list(as.character(1:ntip),
                              as.character(1:ntip))
          diag(va) <- 0
          edgemat <- matrix(ncol=2,nrow=0)
          ## termlens <- diag(v)-colSums(va)
          edgelens <- numeric(0)
          ## maxnode <- ntip
          curnode <- 2*ntip ## one greater than total number of nodes
          ## can we do this in a different order?
          while (nrow(v)>1) {
            mva <- max(va)  ## find pair with max shared evolution
            nextpr <- if (nrow(v)==2) c(1,2) else which(va==mva,arr.ind=TRUE)[1,]
            ## maxnode <- maxnode+1  ## new node
            curnode <- curnode-1
            ## points to both of current identified nodes
            ##   (indexed by names)
            edgemat <- rbind(edgemat,
                             c(curnode,as.numeric(rownames(v)[nextpr[1]])),
                             c(curnode,as.numeric(rownames(v)[nextpr[2]])))
            ## descending edges are amount of *unshared* evolution
            edgelens <- c(edgelens,
                          diag(v)[nextpr]-mva)
            ## this clade has total evolution = shared evolution
            diag(v)[nextpr] <- mva
            ## assign new node name
            rownames(v)[nextpr[1]] <- colnames(v)[nextpr[1]] <- curnode
            ## drop rows/cols from matrix
            v <- v[-nextpr[2],-nextpr[2],drop=FALSE]
            va <- va[-nextpr[2],-nextpr[2],drop=FALSE]
          }
          ## switch order of node numbers to put root in the right place:
          ##  much plotting code seems to assume root = node # (ntips+1)
          ## browser()
          reorder <- FALSE
          if (reorder) {
            nn <- nrow(edgemat)
            nnode <- nn-ntip+1
            newedge <- edgemat
            for (i in 2:nnode) {
              newedge[edgemat==(ntip+i)] <- nn-i+2
            }
            edgemat <- newedge
          }
          list(edgemat=edgemat,
               edgelens=edgelens)
        }
        temptree <- matrix2tree(from)
        ## browser()
        ## add explicit root
        rootnode <- which(tabulate(temptree$edgemat[,2])==0)
        ## add root node to edge matrix and branch lengths
        temptree$edgemat <- rbind(temptree$edgemat, c(0, rootnode))
        temptree$edgelens <- c(temptree$edgelens,NA)
        reorder(phylo4(temptree$edgemat,edge.length=temptree$edgelens,
               tip.label=rownames(from),
               edge.label=from@edge.label,order="unknown"),
                "preorder")
      })


