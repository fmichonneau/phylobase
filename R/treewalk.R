
## matching node labels with node numbers ...
## e.g.
## 14 tips, 13 int nodes
## N04 = nodeLabels[4]
##   <-> node 18
## x = n-nTips(phy)
## so:     n = x+nTips(phy)

getNode <- function(phy, node, missing=c("warn","OK","fail")) {
    missing <- match.arg(missing)
    if (is.numeric(node) && all(floor(node) == node, na.rm=TRUE)) {
        node <- as.integer(node)
    }

    if (is.character(node)) {
        irval <- match(node, labels(phy, "all"))

    }
    else {
        if (is.integer(node)) {
            irval <- match(as.character(node), names(labels(phy, "all")))
        }
        else stop("Node must be a vector of class \'integer\' or \'character\'.")
    }

    ## node numbers
    rval <- names(labels(phy, "all"))[irval]

    rval[node == 0]   <- NA # root ancestor gets special treatment
    rval[is.na(node)] <- NA # return NA for any NA_character_ inputs
    rval <- as.integer(rval)

    ## node labels
    nmNd <- labels(phy, "all")[irval]

    names(rval) <- nmNd
    names(rval)[rval == 0] <- "0" # root ancestor gets special treatment

    ## deal with nodes that don't match
    if (any(is.na(rval))) {
        missnodes <- node[is.na(rval)]
        msg <- paste("Some nodes are missing from tree: ", paste(missnodes,collapse=", "))
        if (missing=="fail") {
            stop(msg)
        } else if (missing=="warn") {
            warning(msg)
        }
    }
    return(rval)
}


ancestor <- function(phy,node) {
    node2 <- getNode(phy,node)
    ## r <- which(edges(phy)[,2]==node)
    r <- match(node2,edges(phy)[,2])
    return(getNode(phy,edges(phy)[r,1],missing="OK"))
}


children <- function(phy,node) {
    node2 <- getNode(phy,node)
    r <- which(edges(phy)[,1]==node2)
    return(getNode(phy,edges(phy)[r,2]))
}

## get descendants [recursively]
descendants <- function (phy, node, type=c("tips","children","all"))
{
    ## FIXME: allow vector of nodes? (or just let people lapply?)
    type <- match.arg(type)
    if (type=="children") {
        return(children(phy, node))
    }
    node <- getNode(phy, node)
    if (is.na(node)) stop("node ", node, " not found in tree")
    n <- nTips(phy)
    if (node <= n) {
        return(node)
    }

    ## edge matrix must be in preorder for the C function!
    if (phy@order=="preorder") {
        edge <- phy@edge
    } else {
        edge <- reorder(phy, order="preorder")@edge
    }
    ## deal with NA root node
    edge[is.na(edge)] <- 0
    ## extract edge columns
    ancestor <- as.integer(edge[, 1])
    descendant <- as.integer(edge[, 2])
    ## create indicator vector and seed it based on children
    isChild <- rep(0L, length=nrow(edge))
    isChild[ancestor==node] <- 1L
    numEdges <- as.integer(length(isChild))

    ## returned value of isChild will indicate *all* descendants 
    isDescendant <- .C("descendants", isChild, ancestor, descendant,
        numEdges)[[1]]

    l <- descendant[as.logical(isDescendant)]
    if (type=="tips") {
        l <- l[l<=n]
    }

    ## Original pure R implementation of the above
    ## (note that it does not require preorder ordering)
    ##l <- numeric()
    ##d <- children(phy, node)
    ##for (j in d) {
    ##    if (j <= n)
    ##      l <- c(l,j)
    ##    else if (type=="all") l <- c(l,j,
    ##               descendants(phy,j,type="all"))
    ##    else l <- c(l, descendants(phy,j,type=type))
    ##}

    return(getNode(phy, l))
}

siblings <- function(phy, node, include.self=FALSE) {
    v <- children(phy,ancestor(phy,node))
    if (!include.self) v <- v[v!=getNode(phy,node)]
    v
}

## get ancestors (all nodes)
ancestors <- function (phy, node, type=c("all","parent","ALL"))
{
    type <- match.arg(type)
    if (type=="parent") return(ancestor(phy,node))
    oNode <- node <- getNode(phy,node)
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
    if(type == "ALL") res <- c(oNode, res)
    return(getNode(phy,res))
}

MRCA <- function(phy, ...) {
    nodes <- list(...)
    ## if length==1 and first element is a vector,
    ##   use it as the list
    if (length(nodes)==1 && length(nodes[[1]])>1) {
        nodes <- as.list(nodes[[1]])
    }

    ## Correct behavior when the root is part of the nodes
    testNodes <- lapply(nodes, getNode, phy=phy)
    ## BMB: why lapply, not sapply?
    lNodes <- unlist(testNodes)
    if (any(is.na(lNodes)))
      stop("nodes not found in tree: ",paste(names(lNodes)[is.na(lNodes)],
                                             collapse=", "))
    uniqueNodes <- unique(testNodes)
    root <- nTips(phy)+1
    if(root %in% uniqueNodes) {
        res <- getNode(phy, root)
        return(res)
    }
    ## Correct behavior in case of MRCA of identical taxa
    if(length(uniqueNodes) == 1) {
        res <- uniqueNodes[[1]]
        return(res)
    }
    else {
        ancests <- lapply(nodes, ancestors, phy=phy, type="ALL")
        res <- getNode(phy, max(Reduce(intersect, ancests)))
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
    ## FIXME: use extractTree if coming from phylo4d

    ## some checks
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)
    t1 <- getNode(x, node1)
    t2 <- getNode(x, node2)
    if(any(is.na(c(t1,t2)))) stop("wrong node specified")
    if(t1==t2) return(NULL)

    ## main computations
    comAnc <- MRCA(x, t1, t2) # common ancestor
    desComAnc <- descendants(x, comAnc, type="all")
    ancT1 <- ancestors(x, t1, type="all")
    path1 <- intersect(desComAnc, ancT1) # path: common anc -> t1

    ancT2 <- ancestors(x, t2, type="all")
    path2 <- intersect(desComAnc, ancT2) # path: common anc -> t2

    res <- union(path1, path2) # union of the path
    ## add the common ancestor if it differs from t1 or t2
    if(!comAnc %in% c(t1,t2)){
        res <- c(comAnc,res)
    }

    res <- getNode(x, res)

    return(res)
} # end shortestPath



###########
# getEdge
###########
getEdge <- function(phy, node, type=c("descendant", "ancestor"),
    missing=c("warn", "OK", "fail")) {

    if(!identical(class(phy), "phylo4")) phy <- as(phy, "phylo4")

    missing <- match.arg(missing)
    node <- getNode(phy, node, missing)

    type <- match.arg(type)

    ##TODO: should missing arg also apply to tips-as-ancestors case?
    nd <- lapply(node, function(x) {
        if (is.na(x)) {
            res <- NA
        } else {
            res <- switch(type,
                descendant = edgeId(phy)[edges(phy)[,2] %in% x],
                ancestor = edgeId(phy)[edges(phy)[,1] %in% x])
            ## hack to return NA for tip nodes when type='ancestor'
            if(length(res)==0) res <- NA
            names(res) <- rep(x, length(res))
        }   
        names(res) <- rep(x, length(res))
        res
    })  

    return(unlist(unname(nd)))

}
