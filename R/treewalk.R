
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
    v <- children(phy,ancestor(phy,node))
    if (!include.self) v <- v[v!=getnodes(phy,node)]
    v
}

## get ancestors (all nodes)
ancestors <- function (phy, node, which=c("all","parent","ALL"))
{
    which <- match.arg(which)
    if (which=="parent") return(ancestor(phy,node))
    oNode <- node <- getnodes(phy,node)
    if (is.na(node)) stop("node ",node," not found in tree")
    res <- numeric(0)
    n <- nTips(phy)

    ## correct behavior when node==root
    if(node == rootNode(phy)) return(NULL)

    repeat {
        anc <- ancestor(phy,node)
        res <- c(res,anc)
        node <- anc
        if (anc==n+1) break
    }
    if(which == "ALL") res <- c(oNode, res)
    return(getnodes(phy,res))
}

MRCA <- function(phy, ...) {
    nodes <- list(...)
    ## if length==1 and first element is a vector,
    ##   use it as the list
    if (length(nodes)==1 && length(nodes[[1]])>1) {
        nodes <- as.list(nodes[[1]])
    }

    ## Correct behavior when the root is part of the nodes
    testNodes <- lapply(nodes, getnodes, phy=phy)
    ## BMB: why lapply, not sapply?
    lNodes <- unlist(testNodes)
    if (any(is.na(lNodes)))
      stop("nodes not found in tree: ",paste(names(lNodes)[is.na(lNodes)],
                                             collapse=", "))
    uniqueNodes <- unique(testNodes)
    root <- nTips(phy)+1
    if(root %in% uniqueNodes) {
        res <- getnodes(phy, root)
        return(res)
    }
    ## Correct behavior in case of MRCA of identical taxa
    if(length(uniqueNodes) == 1) {
        res <- uniqueNodes[[1]]
        return(res)
    }
    else {
        ancests <- lapply(nodes, ancestors, phy=phy, which="ALL")
        res <- getnodes(phy, max(Reduce(intersect, ancests)))
        return(res)
    }
} # end MRCA





###############
# shortestPath
###############
shortestPath <- function(phy, node1, node2){
    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(phy, "phylo4")

    ## come checks
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)
    t1 <- getnodes(x, node1)
    t2 <- getnodes(x, node2)
    if(any(is.na(c(t1,t2)))) stop("wrong node specified")
    if(t1==t2) return(NULL)

    ## main computations
    comAnc <- MRCA(x, t1, t2) # common ancestor
    desComAnc <- descendants(x, comAnc, which="all")
    ancT1 <- ancestors(x, t1, which="all")
    path1 <- intersect(desComAnc, ancT1) # path: common anc -> t1

    ancT2 <- ancestors(x, t2, which="all")
    path2 <- intersect(desComAnc, ancT2) # path: common anc -> t2

    res <- union(path1, path2) # union of the path
    ## add the common ancestor if it differs from t1 or t2
    if(!comAnc %in% c(t1,t2)){
        res <- c(comAnc,res)
    }

    res <- getnodes(x, res)

    return(res)
} # end shortestPath





###########
# getedges
###########
getedges <- function(phy, node){

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(phy, "phylo4")

    ## come checks
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)
    node <- getnodes(x, node)
    if(any(is.na(node))) stop("wrong node specified")
    root <- getnodes(x, nTips(x)+1)
    node[node==root] <- NA

    ## main computations
    E <- x@edge
    res <- match(node, E[,2])
    names(res) <- names(node)

    return(res)
} # end getedges
