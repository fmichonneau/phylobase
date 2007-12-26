
## matching node labels with node numbers ...
## e.g.
## 14 tips, 13 int nodes
## N04 = NodeLabels[4]
##   <-> node 18
## x = n-nTips(phy)
## so:     n = x+nTips(phy)

getNodeByLabel <- function(phy,lab) {
    nt <- phylobase::nTips(phy)
    tipmatch <- match(lab,labels(phy))
    ifelse(!is.na(tipmatch),
           tipmatch,
           if (!hasNodeLabels(phy)) NA else {
               nt+match(lab,NodeLabels(phy))
           })
}
       
getLabelByNode <- function(phy,num) {
    nt <- phylobase::nTips(phy)
    ifelse(num<=nt,
           labels(phy)[num],
           ifelse(num<=nt+nNodes(phy)-1,
                  if (!hasNodeLabels(phy)) NA else {
                      ## pmax to avoid error from negative indices
                      NodeLabels(phy)[pmax(0,num-nt)]
                  }))
}

getAncest <- function(phy,node) {
    if (is.character(node)) node <- getNodeByLabel(phy,node)
    r <- which(phy@edge[,2]==node)
    return(phy@edge[r,1])
}


getDescend <- function(phy,node) {
    if (is.character(node)) node <- getNodeByLabel(phy,node)
    r <- which(phy@edge[,1]==node)
    return(phy@edge[r,2])
}

## get descendants (in this version, tips only)
## recursive
allDescend <- function (phy, node) 
{
    if (is.character(node)) node <- getNodeByLabel(phy,node)
    n <- phylobase::nTips(phy)
    if (node <= n) return(labels(phy)[node])
    l <- character()
    d <- getDescend(phy, node)
    for (j in d) {
        if (j <= n) 
          l <- c(l, labels(phy)[as.numeric(j)])
        else l <- c(l, allDescend(phy, j))
    }
    return(l)
}

## get ancestors (all nodes)
allAncest <- function (phy, node) 
{
    if (is.character(node)) node <- getNodeByLabel(phy,node)
    res <- numeric(0)
    n <- phylobase::nTips(phy)
    repeat {
        anc <- getAncest(phy,node)
        res <- c(res,anc)
        node <- anc
        if (anc==n+1) break
    }
    return(res)
}

MRCA <- function(phy, ...) {
    nodes <- list(...)
    ## if length==1 and first element is a vector,
    ##   use it as the list
    if (length(nodes)==1 && length(nodes[[1]])>1) {
        nodes <- as.list(nodes[[1]])
    }
    ancests <- lapply(nodes,allAncest,phy=phy)
    max(Reduce(intersect,ancests))
}
