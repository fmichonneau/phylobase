## matching node labels with node numbers ...
## e.g.
## 14 tips, 13 int nodes
## N04 = nodeLabels[4]
##   <-> node 18
## x = n-nTips(phy)
## so:     n = x+nTips(phy)



#' node and edge look-up functions
#' 
#' Functions for retrieving node and edge IDs (possibly with corresponding
#' labels) from a phylogenetic tree.
#' 
#' \code{getNode} and \code{getEdge} are primarily intended for looking up the
#' IDs either of nodes themselves or of edges associated with those nodes. Note
#' that they behave quite differently. With \code{getNode}, any input nodes are
#' looked up against tree nodes of the specified type, and those that match are
#' returned as numeric node IDs with node labels (if they exist) as element
#' names. With \code{getEdge}, any input nodes are looked up against edge ends
#' of the specified type, and those that match are returned as character edge
#' IDs with the corresponding node ID as element names.
#' 
#' If \code{missing} is \dQuote{warn} or \dQuote{OK}, \code{NA} is returned for
#' any nodes that are unmatched for the specified type. This can provide a
#' mechanism for filtering a set of nodes or edges.
#' 
#' \code{nodeId} provides similar output to \code{getNode} in the case when no
#' node is supplied, but it is faster and returns an unnamed vector of the
#' numeric IDs of all nodes of the specified node type.  Similarly,
#' \code{edgeId} simply returns an unnamed vector of the character IDs of all
#' edges for which the descendant node is of the specified node type.
#' 
#' @aliases getNode getEdge nodeId nodeId,phylo4-method edgeId
#' edgeId,phylo4-method
#' @param x a \linkS4class{phylo4} object (or one inheriting from
#' \linkS4class{phylo4}, e.g. a \linkS4class{phylo4d} object)
#' @param node either an integer vector corresponding to node ID numbers, or a
#' character vector corresponding to node labels; if missing, all nodes
#' appropriate to the specified type will be returned by \code{getNode}, and
#' all edges appropriate to the specified type will be returned by
#' \code{getEdge}.
#' @param type (\code{getNode}) specify whether to return nodes matching "all"
#' tree nodes (default), only "tip" nodes, or only "internal" nodes;
#' (\code{nodeId, edgeId}) specify whether to return "all" tree nodes, or only
#' those corresponding to "tip", "internal", or "root" nodes; (\code{getEdge})
#' specify whether to look up edges based on their descendant node
#' ("descendant") or ancestral node ("ancestor")
#' @param missing what to do if some requested node IDs or names are not in the
#' tree: warn, do nothing, or stop with an error
#' @return \item{list("getNode")}{returns a named integer vector of node IDs,
#' in the order of input nodes if provided, otherwise in nodeId order}
#' \item{list("getEdge")}{returns a named character vector of edge IDs, in the
#' order of input nodes if provide, otherwise in nodeId order}
#' \item{list("nodeId")}{returns an unnamed integer vector of node IDs, in
#' ascending order} \item{list("getEdge")}{returns an unnamed character vector
#' of edge IDs, in edge matrix order}
#' @keywords misc
#' @examples
#' 
#'   data(geospiza)
#'   nodeLabels(geospiza) <- LETTERS[1:nNodes(geospiza)]
#'   plot(as(geospiza, "phylo4"), show.node.label=TRUE)
#'   getNode(geospiza, 18)
#'   getNode(geospiza, "D")
#'   getEdge(geospiza, "D")
#'   getEdge(geospiza, "D", type="ancestor")
#' 
#'   ## match nodes only to tip nodes, flagging invalid cases as NA
#'   getNode(geospiza, c(1, 18, 999), type="tip", missing="OK")
#' 
#'   ## get all edges that descend from internal nodes
#'   getEdge(geospiza, type="ancestor")
#' 
#'   ## identify an edge from its terminal node
#'   getEdge(geospiza, c("olivacea", "B", "fortis"))
#'   getNode(geospiza, c("olivacea", "B", "fortis"))
#'   edges(geospiza)[c(26, 1, 11),]
#' 
#'   ## quickly get all tip node IDs and tip edge IDs
#'   nodeId(geospiza, "tip")
#'   edgeId(geospiza, "tip")
#' 
getNode <- function(x, node, type=c("all", "tip", "internal"),
    missing=c("warn","OK","fail")) {

    type <- match.arg(type)
    missing <- match.arg(missing)

    ## if missing node arg, get all nodes of specified type
    if (missing(node)) {
        node <- nodeId(x, type)
    }

    if (length(node) == 0) {
      rval <- integer(0)
      names(rval) <- character(0)
      return(rval)
    }
    
    lblTmp <- labels(x, type)
    
    ## match node to tree
    if (is.character(node)) {
        ndTmp <- paste("^\\Q", node, "\\E$", sep="")        
        irval <- lapply(ndTmp, function(ND) {
            grep(ND, lblTmp, perl=TRUE)
        })
        irvalL <- sapply(irval, length)
        irval[irvalL == 0] <- 0
        irval <- unlist(irval)
    } else if (is.numeric(node) && all(floor(node) == node, na.rm=TRUE)) {
        irval <- match(as.character(node), names(lblTmp))
    } else {
        stop("Node must be a vector of class \'integer\' or \'character\'.")
    }

    ## node numbers
    rval <- names(lblTmp)[irval]
    rval[is.na(node)] <- NA # return NA for any NA_character_ inputs, not needed but ensure rval has correct length
    rval <- as.integer(rval)

    ## node labels
    nmNd <- lblTmp[irval]
    names(rval) <- nmNd
    
    ## deal with nodes that don't match
    if (any(is.na(rval))) {
        missnodes <- node[is.na(rval)]
        msg <- paste("Some nodes not found among", type, "nodes in tree:",
            paste(missnodes,collapse=", "))
        if (missing=="fail") {
            stop(msg)
        } else if (missing=="warn") {
            warning(msg)
        }
    }
    return(rval)
}


#' tree traversal and utility functions
#' 
#' Functions for describing relationships among phylogenetic nodes (i.e.
#' internal nodes or tips).
#' 
#' \code{ancestors} and \code{descendants} can take \code{node} vectors of
#' arbitrary length, returning a list of output vectors if the number of valid
#' input nodes is greater than one. List element names are taken directly from
#' the input node vector.
#' 
#' If any supplied nodes are not found in the tree, the behavior currently
#' varies across functions. Invalid nodes are automatically omitted by
#' \code{ancestors} and \code{descendants}, with a warning.  \code{ancestor}
#' will return \code{NA} for any invalid nodes, with a warning. Both
#' \code{children} and \code{siblings} will return an empty vector, again with
#' a warning. In contrast, \code{MRCA} and \code{shortestPath} will throw an
#' immediate error if any input nodes are invalid.
#' 
#' @aliases children descendants ancestor ancestors siblings MRCA shortestPath
#' sumEdgeLength sumEdgeLength,phylo4-method
#' @param phy a \linkS4class{phylo4} object (or one inheriting from
#' \linkS4class{phylo4}, e.g. a \linkS4class{phylo4d} object)
#' @param x a \linkS4class{phylo4} object (or one inheriting from
#' \linkS4class{phylo4}, e.g. a \linkS4class{phylo4d} object)
#' @param node either an integer corresponding to a node ID number, or a
#' character corresponding to a node label; for \code{ancestors} and
#' \code{descendants}, this may be a vector of multiple node numbers or names
#' @param type (\code{ancestors}) specify whether to return just direct
#' ancestor ("parent"), all ancestor nodes ("all"), or all ancestor nodes
#' including self ("ALL"); (\code{descendants}) specify whether to return just
#' direct descendants ("children"), all extant descendants ("tips"), or all
#' descendant nodes ("all")
#' @param include.self whether to include self in list of siblings
#' @param \dots a list of node numbers or names, or a vector of node numbers or
#' names
#' @param node1 a node number (or name)
#' @param node2 a node number (or name)
#' @return \item{list("ancestors")}{ return a named vector (or a list of such
#' vectors in the case of multiple input nodes) of the ancestors and
#' descendants of a node}\item{ and }{ return a named vector (or a list of such
#' vectors in the case of multiple input nodes) of the ancestors and
#' descendants of a node}\item{list("descendants")}{ return a named vector (or
#' a list of such vectors in the case of multiple input nodes) of the ancestors
#' and descendants of a node} \item{list("ancestor")}{ \code{ancestor} is
#' analogous to \code{ancestors(\dots{}, type="parent")} (i.e. direct ancestor
#' only), but returns a single concatenated vector in the case of multiple
#' input nodes; \code{children} is analogous to \code{descendants(\dots{},
#' type="children")} (i.e. direct descendants only), but is not currently
#' intended to be used with multiple input nodes }\item{ and }{ \code{ancestor}
#' is analogous to \code{ancestors(\dots{}, type="parent")} (i.e. direct
#' ancestor only), but returns a single concatenated vector in the case of
#' multiple input nodes; \code{children} is analogous to
#' \code{descendants(\dots{}, type="children")} (i.e. direct descendants only),
#' but is not currently intended to be used with multiple input nodes
#' }\item{list("children")}{ \code{ancestor} is analogous to
#' \code{ancestors(\dots{}, type="parent")} (i.e. direct ancestor only), but
#' returns a single concatenated vector in the case of multiple input nodes;
#' \code{children} is analogous to \code{descendants(\dots{}, type="children")}
#' (i.e. direct descendants only), but is not currently intended to be used
#' with multiple input nodes } \item{list("siblings")}{ returns sibling nodes
#' (children of the same parent)} \item{list("MRCA")}{ returns the most recent
#' common ancestor of two or more nodes} \item{list("shortestPath")}{ returns
#' the nodes of the shortest path from one node to another (excluding
#' \code{node1} and \code{node2})} \item{list("sumEdgeLength")}{ returns the
#' sum of branch length for branches starting at nodes provided}
#' @note \code{MRCA} is uppercase to avoid conflict with \code{mrca} in ape
#' @seealso \code{\link[ape]{mrca}}, in the ape package, gives a list of all
#' subtrees
#' @keywords misc
#' @examples
#' 
#'   data(geospiza)
#'   nodeLabels(geospiza) <- LETTERS[1:nNodes(geospiza)]
#'   plot(as(geospiza, "phylo4"), show.node.label=TRUE)
#'   ancestor(geospiza, "E")
#'   children(geospiza, "C")
#'   descendants(geospiza, "D", type="tips")
#'   descendants(geospiza, "D", type="all")
#'   ancestors(geospiza, "D")
#'   MRCA(geospiza, "conirostris", "difficilis", "fuliginosa")
#'   MRCA(geospiza, "olivacea", "conirostris")
#' 
#'   ## shortest path between 2 nodes
#'   shortestPath(geospiza, "fortis", "fuliginosa")
#'   shortestPath(geospiza, "F", "L")
#' 
#'   ## branch length from a tip to the root
#'   sumEdgeLength(geospiza, ancestors(geospiza, "fortis", type="ALL"))
ancestor <- function(phy,node) {
    node2 <- getNode(phy,node)
    ## r <- which(edges(phy)[,2]==node)
    r <- match(node2,edges(phy)[,2])
    return(getNode(phy,edges(phy)[r,1],missing="OK"))
}


children <- function(phy,node) {
    node2 <- getNode(phy,node)
    r <- which(edges(phy)[,1]==node2)
    getNode(phy,edges(phy)[r,2])
}

## get descendants [recursively]
descendants <- function (phy, node, type=c("tips","children","all")) {
    type <- match.arg(type)

    ## look up nodes, warning about and excluding invalid nodes
    oNode <- node
    node <- getNode(phy, node, missing="warn")
    isValid <- !is.na(node)
    node <- as.integer(node[isValid])

    if (type == "children") {
        res <- lapply(node, function(x) children(phy, x))
        ## if just a single node, return as a single vector
        if (length(res)==1) res <- res[[1]]
    } else {
        ## edge matrix must be in preorder for the C function!
        if (phy@order=="preorder") {
            edge <- phy@edge
        } else {
            edge <- reorder(phy, order="preorder")@edge
        }
        ## extract edge columns
        ancestor <- as.integer(edge[, 1])
        descendant <- as.integer(edge[, 2])

        ## return indicator matrix of ALL descendants (including self)
        isDes <- .Call("descendants", node, ancestor, descendant)
        storage.mode(isDes) <- "logical"

        ## for internal nodes only, drop self (not sure why this rule?)
        int.node <- intersect(node, nodeId(phy, "internal"))
        isDes[cbind(match(int.node, descendant),
            match(int.node, node))] <- FALSE
        ## if only tips desired, drop internal nodes
        if (type=="tips") {
            isDes[descendant %in% nodeId(phy, "internal"),] <- FALSE
        }
        ## res <- lapply(seq_along(node), function(n) getNode(phy,
        ##     descendant[isDes[,n]]))
        res <- getNode(phy, descendant[isDes[, seq_along(node)]])
    }
    ## names(res) <- as.character(oNode[isValid])

    res

    ## Original pure R implementation of the above
    ## (note that it does not require preorder ordering)
    ##n <- nTips(phy)
    ##if (node <= n) {
    ##    return(node)
    ##}
    ##l <- numeric()
    ##d <- children(phy, node)
    ##for (j in d) {
    ##    if (j <= n)
    ##      l <- c(l,j)
    ##    else if (type=="all") l <- c(l,j,
    ##               descendants(phy,j,type="all"))
    ##    else l <- c(l, descendants(phy,j,type=type))
    ##}
}

siblings <- function(phy, node, include.self=FALSE) {
    v <- children(phy,ancestor(phy,node))
    if (!include.self) v <- v[v!=getNode(phy,node)]
    v
}

## get ancestors (all nodes)
ancestors <- function (phy, node, type=c("all","parent","ALL")) {

    type <- match.arg(type)

    ## look up nodes, warning about and excluding invalid nodes
    oNode <- node
    node <- getNode(phy, node, missing="warn")
    isValid <- !is.na(node)
    node <- as.integer(node[isValid])

    if (length(node) == 0) {
      return(NA)
    }
    
    if (type == "parent") {
        res <- lapply(node, function(x) ancestor(phy, x))
    } else {
        ## edge matrix must be in postorder for the C function!
        if (phy@order=="postorder") {
            edge <- phy@edge
        } else {
            edge <- reorder(phy, order="postorder")@edge
        }
        ## extract edge columns
        ancestor <- as.integer(edge[, 1])
        descendant <- as.integer(edge[, 2])

        ## return indicator matrix of ALL ancestors (including self)
        isAnc <- .Call("ancestors", node, ancestor, descendant)
        storage.mode(isAnc) <- "logical"

        ## drop self if needed
        if (type=="all") {
            isAnc[cbind(match(node, descendant), seq_along(node))] <- FALSE
        }
        res <- lapply(seq_along(node), function(n) getNode(phy,
            descendant[isAnc[,n]]))
    }
    names(res) <- as.character(oNode[isValid])

    ## if just a single node, return as a single vector
    if (length(res)==1) res <- res[[1]]
    res

    ## Original pure R implementation of the above
    ## (note that it does not require preorder ordering)
    ##if (node == rootNode(phy))
    ##    return(NULL)
    ##repeat {
    ##    anc <- ancestor(phy, node)
    ##    res <- c(res, anc)
    ##    node <- anc
    ##    if (anc == n + 1)
    ##        break
    ##}
}

MRCA <- function(phy, ...) {
    nodes <- list(...)
    ## if length==1 and first element is a vector,
    ##   use it as the list
    if (length(nodes)==1 && length(nodes[[1]])>1) {
        nodes <- as.list(nodes[[1]])
    }

    ## Correct behavior when the root is part of the nodes
    testNodes <- lapply(nodes, getNode, x=phy)
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

  ## conversion from phylo, phylo4 and phylo4d
  if (class(phy) == "phylo4d") {
    x <- extractTree(phy)
  }
  else if (class(phy) != "phylo4"){
    x <- as(phy, "phylo4")
  }

    ## some checks
    ## if (is.character(checkval <- checkPhylo4(x))) stop(checkval) # no need
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
getEdge <- function(x, node, type=c("descendant", "ancestor"),
    missing=c("warn", "OK", "fail")) {

    if(!identical(class(x), "phylo4")) x <- as(x, "phylo4")

    type <- match.arg(type)
    missing <- match.arg(missing)
    if (missing(node)) {
        if (type=="descendant") {
            node <- nodeId(x, "all")
        } else if (type=="ancestor") {
            node <- nodeId(x, "internal")
        }
    }

    node.id <- getNode(x, node, missing="OK")

    nd <- lapply(node.id, function(nid) {
        if (is.na(nid)) {
            res <- NA
        } else {
            res <- switch(type,
                descendant = edgeId(x)[edges(x)[,2] %in% nid],
                ancestor = edgeId(x)[edges(x)[,1] %in% nid])
            ## hack to return NA for tip nodes when type='ancestor'
            if(length(res)==0) res <- NA
            names(res) <- rep(nid, length(res))
        }
        names(res) <- rep(nid, length(res))
        res
    })

    ## warn or stop if necessary
    is.missing <- is.na(nd)
    if (missing!="OK" && any(is.missing)) {
        msg <- paste("Not all nodes are ", type, "s in this tree: ",
            paste(node[is.missing], collapse=", "), sep="")
        if (missing=="fail") {
            stop(msg)
        } else if (missing=="warn") {
            warning(msg)
        }
    }

    return(unlist(unname(nd)))

}
