getNodeByLabel <- function(phy,lab) {
    tipmatch <- match(lab,labels(phy))
    if (!is.na(tipmatch)) return(tipmatch)
    if (!hasNodeLabels(phy)) return(NULL)
    which(NodeLabels(phy)==lab)+phylobase::nTips(phy)
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
