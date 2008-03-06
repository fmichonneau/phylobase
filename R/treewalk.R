
## matching node labels with node numbers ...
## e.g.
## 14 tips, 13 int nodes
## N04 = NodeLabels[4]
##   <-> node 18
## x = n-nTips(phy)
## so:     n = x+nTips(phy)

getnodes <- function(phy,x) {
    if (is.numeric(x) && all(floor(x)==x)) {
        x <- as.integer(x)
    }
    if (is.character(x)) {
        ## old getNodeByLabel()
        nt <- nTips(phy)
        tipmatch <- match(x,labels(phy,"all"))
        vals <- ifelse(!is.na(tipmatch),
                       tipmatch,
                       if (!hasNodeLabels(phy)) { NA } else {
                           nt+match(x,NodeLabels(phy))
                       })
        names(vals) <- x
        return(vals)
    } else if (is.integer(x)) {
        ## old getLabelByNode
        nt <- nTips(phy)
        vals <- ifelse(x<=nt,  ## tips
                       labels(phy,"all")[x], 
                       ifelse(x<=nt+nNodes(phy),
                              if (!hasNodeLabels(phy)) { NA }
                              else {
                                  NodeLabels(phy)[pmax(0,x-nt)]
                              },NA))
        ## pmax above to avoid error from negative indices
        names(x) <- vals
        return(x)
    } else stop("x must be integer or character")
}


parent <- function(phy,node) {
    node <- getnodes(phy,node)
    r <- which(phy@edge[,2]==node)
    return(getnodes(phy,phy@edge[r,1]))
}


children <- function(phy,node) {
    node <- getnodes(phy,node)
    r <- which(phy@edge[,1]==node)
    return(getnodes(phy,phy@edge[r,2]))
}

## get descendants [recursively]
descendants <- function (phy, node, which=c("tips","all"))
{
    ## FIXME: allow vector of nodes?
    which <- match.arg(which)
    if (which=="all") stop('which="all" not yet working')
    node <- getnodes(phy,node)
    if (is.na(node)) stop("node ",node," not found in tree")
    n <- nTips(phy)
    if (node <= n) return(labels(phy,"all")[node])
    l <- numeric()
    d <- children(phy, node)
    for (j in d) {
        if (j <= n)
          l <- c(l,j)
        else if (which=="all") l <- c(l,node,descendants(phy,j,which="all"))
        else l <- c(l, descendants(phy,j))
##           l <- c(l, labels(phy,"tips")[as.numeric(j)])
##         else if (which=="all") l <- c(l, names(node),
##                    names(descendants(phy, j, which="all")))
##         else l <- c(l, names(descendants(phy, j)))
    }
    return(getnodes(phy,l))
}

## get ancestors (all nodes)
ancestors <- function (phy, node) 
{
    node <- getnodes(phy,node)
    if (is.na(node)) stop("node ",node," not found in tree")
    res <- numeric(0)
    n <- nTips(phy)
    repeat {
        anc <- parent(phy,node)
        res <- c(res,anc)
        node <- anc
        if (anc==n+1) break
    }
    return(getnodes(phy,res))
}

MRCA <- function(phy, ...) {
    nodes <- list(...)
    ## if length==1 and first element is a vector,
    ##   use it as the list
    if (length(nodes)==1 && length(nodes[[1]])>1) {
        nodes <- as.list(nodes[[1]])
    }
    ancests <- lapply(nodes,ancestors,phy=phy)
    getnodes(phy,max(Reduce(intersect,ancests)))
}
