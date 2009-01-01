
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
        temptree$edgemat <- rbind(temptree$edgemat,c(NA,rootnode))
        temptree$edgelens <- c(temptree$edgelens,NA)
        phylo4(temptree$edgemat,edge.length=temptree$edgelens,
               tip.label=rownames(from),
               edge.label=from@edge.label,order=from@order)
      })


