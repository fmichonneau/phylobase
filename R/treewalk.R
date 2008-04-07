
## matching node labels with node numbers ...
## e.g.
## 14 tips, 13 int nodes
## N04 = nodeLabels[4]
##   <-> node 18
## x = n-nTips(phy)
## so:     n = x+nTips(phy)

getnodes <- function(phy,node) {
    if (is.numeric(node) && all(floor(node)==node)) {
        node <- as.integer(node)
    }
    if (is.character(node)) {
        ## old getNodeByLabel()
        nt <- nTips(phy)
        tipmatch <- match(node,labels(phy,"all"))
        vals <- ifelse(!is.na(tipmatch),
                       tipmatch,
                       if (!hasNodeLabels(phy)) { NA } else {
                           nt+match(node,nodeLabels(phy))
                       })
        names(vals) <- node
        return(vals)
    } else if (is.integer(node)) {
        ## old getLabelByNode
        nt <- nTips(phy)
        vals <- ifelse(node<=nt,  ## tips
                       labels(phy,"all")[node], 
                       ifelse(node<=nt+nNodes(phy),
                              if (!hasNodeLabels(phy)) { NA }
                              else {
                                  nodeLabels(phy)[pmax(0,node-nt)]
                              },NA))
        ## pmax above to avoid error from negative indices
        names(node) <- vals
        return(node)
    } else stop("node must be integer or character")
}


ancestor <- function(phy,node) {
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
descendants <- function (phy, node, which=c("tips","children","all"))
{
    ## FIXME: allow vector of nodes? (or just let people lapply?)
    which <- match.arg(which)
    if (which=="children") return(children(phy,node))
    node <- getnodes(phy,node)
    if (is.na(node)) stop("node ",node," not found in tree")
    n <- nTips(phy)
    if (node <= n) return(labels(phy,"all")[node])
    l <- numeric()
    d <- children(phy, node)
    for (j in d) {
        if (j <= n)
          l <- c(l,j)
        else if (which=="all") l <- c(l,j,
                   descendants(phy,j,which="all"))
        else l <- c(l, descendants(phy,j,which=which))
    }
    return(getnodes(phy,l))
}

siblings <- function(phy, node, include.self=FALSE) {
    v <- children(phy,parent(phy,node))
    if (!include.self) v <- v[v!=getnodes(phy,node)]
    v
}
    
## get ancestors (all nodes)
ancestors <- function (phy, node, which=c("all","parent")) 
{
    which <- match.arg(which)
    if (which=="parent") return(ancestor(phy,node))
    node <- getnodes(phy,node)
    if (is.na(node)) stop("node ",node," not found in tree")
    res <- numeric(0)
    n <- nTips(phy)
    repeat {
        anc <- ancestor(phy,node)
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
