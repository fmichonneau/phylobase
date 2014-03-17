### This file contains the methods and accessors for phylo4(d) objects
### The file is organized in sections:

### 1. Tip accessors
###  1.1. nTips()
###  1.2. depthTips()

### 2. Node accessors
###  2.1. nNodes()
###  2.2. nodeType()
###  2.3. nodeId()
###  2.4. nodeDepth()

### 3. Edge accessors
###  3.1. nEdges()
###  3.2. edges()
###  3.3. edgeOrder()
###  3.4. edgeId()
###  3.5. hasEdgeLength()
###  3.6. edgeLength()
###  3.7. edgeLength() <-
###  3.8. sumEdgeLength()

### 4. Root accessors
###  4.1. isRooted()
###  4.2. rootNode()
###  4.3. rootNode() <-

### 5. Label accessors
###  5.1. labels()
###  5.2. labels() <-
###  5.3. hasNodeLabels()
###  5.4. nodeLabels()
###  5.5. nodeLabels() <-
###  5.6. tipLabels()
###  5.7. tipLabels() <-
###  5.8. hasEdgeLabels()
###  5.9. edgeLabels()
###  5.10. edgeLabels() <-

### 6. Displaying functions
###  6.1. printphylo4()
###  6.2. print()
###  6.3  show()
###  6.4. names()
###  6.5. head()
###  6.6. tail()
###  6.7. summary()

### 7. Ordering
###  7.1. orderIndex()
###  7.2. reorder()

### 8. Tree properties
###  8.1. isUltrametric()

#' 
#' @name phylo4-accessors
#' @aliases nNodes nNodes-methods nNodes,phylo4-method nTips nTips-methods
#' nTips,phylo4-method nTips,phylo-method depthTips depthTips-methods
#' depthTips,phylo4-method depthTips,phylo4d-method edges edges-methods
#' edges,phylo4-method nEdges nEdges-methods nEdges,phylo4-method nodeDepth
#' nodeDepth-methods nodeDepth,phylo4-method edgeOrder edgeOrder,phylo4-method
#' hasEdgeLength hasEdgeLength-methods hasEdgeLength,phylo4-method edgeLength
#' edgeLength-methods edgeLength,phylo4-method edgeLength<-
#' edgeLength<-,phylo4-method edgeLength<-,phylo4,ANY-method nodeType
#' nodeType,phylo4-method isRooted isRooted-methods isRooted,phylo4-method
#' rootEdge rootEdge-methods rootEdge,phylo4-method rootNode rootNode-methods
#' rootNode,phylo4-method rootNode<- rootNode<-,phylo4-method isUltrametric
#' isUltrametric-methods isUltrametric,phylo4-method
#' @docType methods
#' @param x a phylo4/phylo4d object
#' @param node which edge to extract (indexed by descendant node)
#' @param value a vector of edge lengths or a node number
#' @param use.names Should the names of \code{value} be used to match edge
#' lengths provided?
#' @param drop.root logical: drop root row from edge matrix?
#' @param tol tolerance in rounding error to determine whether the tree is
#' ultrametric
#' @param \dots additional parameters passed (currently ignored)
#' @section Methods: \describe{ \item{nTips}{\code{signature(object="phylo4")}:
#' number of tips}
#' 
#' \item{depthTips}{\code{signature(object="phylo4")}: distance between the
#' tips and the root}
#' 
#' \item{nNodes}{\code{signature(object="phylo4")}: number of internal nodes}
#' 
#' \item{nEdges}{\code{signature(object = "phylo4")}: number of edges}
#' 
#' \item{edges}{\code{signature(object = "phylo4")}: returns the edge matrix}
#' 
#' \item{edgeOrder}{\code{signature(object = "phylo4")}: returns the order in
#' which the edges are stored}
#' 
#' \item{hasEdgeLength}{\code{signature(object = "phylo4")}: whether tree has
#' edge (branch) lengths}
#' 
#' \item{edgeLength}{\code{signature(object = "phylo4")}: edge (branch) lengths
#' (or NAs if missing) ordered according to the edge matrix}
#' 
#' \item{nodeType}{\code{signature(object = "phylo4")}: named vector which has
#' the type of node (internal, tip, root) for value, and the node number for
#' name}
#' 
#' \item{nodeDepth}{\code{signature(object = "phylo4")}: named vector which
#' gives the distance between nodes and the root}
#' 
#' \item{isRooted}{\code{signature(object = "phylo4")}: whether tree is rooted
#' (i.e. has explicit root edge defined \emph{or} root node has <= 2
#' descendants)}
#' 
#' \item{rootEdge}{\code{signature(object = "phylo4")}: root edge}
#' 
#' \item{isUltrametric}{\code{signature(object = "phylo4")}: whether the tree
#' is ultrametric} }
#' @keywords methods
#' @examples
#' 
#' data(geospiza)
#' edgeLength(geospiza, 5)
#' edgeLength(geospiza, "olivacea")
#' edgeLength(geospiza, 5:7)
#' 

#########################################################
### Tip accessors
#########################################################

setMethod("nTips", signature(x="phylo4"), function(x) {
    E <- edges(x)
    if(nrow(E) == 0)
        return(0)
    else {
        ## at this time NAs are not allowed in edge matrix
        ## sum(tabulate(E[, 1]) == 0)
        nTipsFastCpp(E[, 1])
    }        
})

## hack to ensure ape compatibility
setMethod("nTips", signature(x="phylo"),
 function(x) {
     Ntip(x)
})

setMethod("depthTips", signature(x="phylo4"), function(x) {
  nodeDepth(x, 1:nTips(x))
})

#########################################################
### Node accessors
#########################################################

setMethod("nNodes", signature(x="phylo4"), function(x) {
    E <- edges(x, drop.root=TRUE)
    if(nrow(E) == 0) {
        return(0)
    } else {
        return(length(unique(E[, 1])))
    }
})

setMethod("nodeType", signature(x="phylo4"),
 function(x) {
    if(nTips(x) == 0)
        return(NULL)
    else {
        ## strip out the root ancestor
        nodesVect <- as.vector(edges(x))
        nodesVect <- nodesVect[nodesVect != 0]
        ## get a sorted list of the unique nodes
        listNodes <- sort(unique(nodesVect))
        t <- rep("internal", length(listNodes)) # FM: internal is default (I think it's safer)
        names(t) <- listNodes

        ## node number of real internal nodes
        iN <- names(table(edges(x)[,1]))
        ## node number that are not internal nodes (ie that are tips)
        tN <- names(t)[!names(t) %in% iN]
        t[tN] <- "tip"

        ## if the tree is rooted
        if(isRooted(x)) t[rootNode(x)] <- "root"

        return(t)
    }
})

# return node IDs (or a subset thereof) in ascending order
setMethod("nodeId", signature(x="phylo4"),
 function(x, type=c("all",
    "tip","internal","root")) {

     type <- match.arg(type)
     E <- edges(x)

     ## Note: this implementation will still work even if tips are not
     ## 1:nTips and nodes are not (nTips+1):nNodes
     nid <- switch(type,
         ## all nodes appear at least once in the edge matrix
         all = unique(as.vector(E)[as.vector(E) != 0]),
         ## tips are nodes that do not appear in the ancestor column
         tip = setdiff(E[, 2], E[, 1]),
         ## internals are nodes that *do* appear in the ancestor column
         internal = unique(E[E[, 1] != 0, 1]),
         ## roots are nodes that have NA as ancestor
         root = if (!isRooted(x)) NA else unname(E[E[, 1] == 0, 2]))

     return(sort(nid))

})

## nodeId
setGeneric("nodeIdCpp", function(x, type=c("all", "tip", "internal",
    "root")) {
    standardGeneric("nodeIdCpp")
})

setMethod("nodeIdCpp", signature(x="phylo4"),
          function(x, type=c("all", "tip", "internal", "root")) {
              type <- match.arg(type)
              E <- edges(x)
              nid <- switch(type,
                            all = getAllNodesFast(x@edge, isRooted(x)),
                            tip = tipsFast(x@edge[,1]),
                            internal = setdiff(getAllNodesFast(x@edge, isRooted(x)), tipsFast(x@edge[,1])),
                            root = if (!isRooted(x)) NA else unname(E[E[, 1] == 0, 2]))
              nid
          })

setMethod("nodeDepth", signature(x="phylo4"),
  function(x, node) {
    if (!hasEdgeLength(x))
      return(NULL)
    else {
      node <- getNode(x, node, missing="warn")
      node <- node[!is.na(node)]
      res <- sapply(node, function(n)
                    sumEdgeLength(x, ancestors(x, n, "ALL")))
      if (length(res) == 1) {
        res <- res[[1]]
        names(res) <- names(node)
      }      
      res
    }
})

#########################################################
### Edge accessors
#########################################################

setMethod("nEdges", signature(x="phylo4"),
 function(x) {
    nrow(x@edge)
})

# return edge matrix in its current order
setMethod("edges", signature(x="phylo4"),
 function(x, drop.root=FALSE, ...) {
     e <- x@edge
     if (drop.root) e <- e[e[, 1] != 0, ]
     e
})

setMethod("edgeOrder", signature(x="phylo4"),
 function(x, ...) {
    x@order
})

# return edge IDs (or a subset thereof) in edge matrix order
setMethod("edgeId", signature(x="phylo4"),
 function(x, type=c("all", "tip",
    "internal", "root")) {
    type <- match.arg(type)
    edge <- edges(x)
    if (type=="tip") {
        isTip <- !(edge[, 2] %in% edge[, 1])
        edge <- edge[isTip, , drop=FALSE]
    } else if (type=="internal") {
        isInt <- (edge[, 2] %in% edge[, 1])
        edge <- edge[isInt, , drop=FALSE]
    } else if (type=="root") {
        isRoot <- edge[, 1] == 0
        edge <- edge[isRoot, , drop=FALSE]
    } # else just use complete edge matrix if type is "all"
    id <- paste(edge[, 1], edge[, 2], sep="-")
    return(id)
})

setMethod("hasEdgeLength", signature(x="phylo4"),
 function(x) {
    !all(is.na(x@edge.length))
})

# return edge lengths in order by edgeIds (same order as edge matrix)
setMethod("edgeLength", signature(x="phylo4"),
 function(x, node) {
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    if (!missing(node)) {
        id <- getEdge(x, node)
    } else {
        id <- edgeId(x, "all")
    }
    elen <- x@edge.length[match(id, names(x@edge.length))]
    names(elen) <- id
    return(elen)
})

setReplaceMethod("edgeLength", signature(x="phylo4"),
 function(x, use.names=TRUE, ..., value) {
    len <- .createEdge(value, x@edge, type="lengths", use.names)
    ## return empty vector if all values are NA
    if (all(is.na(len))) {
        emptyVec <- numeric()
        attributes(emptyVec) <- list(names=character(0))
        x@edge.length <- emptyVec
    } else {
        x@edge.length <- len
    }
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    x
})

setMethod("sumEdgeLength", signature(x="phylo4"),
 function(x, node) {
    if(!hasEdgeLength(x))
        NULL
    else {
        nd <- getNode(x, node)
        iEdges <- which(x@edge[,2] %in% nd)
        sumEdges <- sum(x@edge.length[iEdges],na.rm=TRUE)
        sumEdges
    }
})

#########################################################
### Root accessors
#########################################################

setMethod("isRooted", signature(x="phylo4"),
 function(x) {
    ## hack to avoid failure on an empty object
    if(nTips(x) == 0) return(FALSE)
    any(edges(x)[, 1] == 0)
})

setMethod("rootNode", signature(x="phylo4"),
 function(x) {
    if (!isRooted(x))
        return(NA)
    unname(edges(x)[which(edges(x)[, 1] == 0), 2])
})

setReplaceMethod("rootNode", signature(x="phylo4"),
 function(x, value) {
    stop("Root node replacement not implemented yet")
})

#########################################################
### Label accessors
#########################################################

#' Labels for phylo4/phylo4d objects
#' 
#' Methods for creating, accessing and updating labels in phylo4/phylo4d
#' objects
#' 
#' 
#' In phylo4/phylo4d objects, tips must have labels (that's why there is no
#' method for hasTipLabels), internal nodes and edges can have labels.
#' 
#' Labels must be provided as a vector of class \code{character}. The length of
#' the vector must match the number of elements they label.
#' 
#' The option \code{use.names} allows the user to match a label to a particular
#' node. In this case, the vector must have names that match the node numbers.
#' 
#' The function \code{labels} is mostly intended to be used internally.
#' 
#' @name phylo4-labels
#' @aliases labels<- labels,phylo4-method
#' labels<-,phylo4,ANY,ANY,character-method
#' labels<-,phylo4d,ANY,ANY,character-method hasDuplicatedLabels
#' hasDuplicatedLabels-methods hasDuplicatedLabels,phylo4-method hasNodeLabels
#' hasNodeLabels-methods hasNodeLabels,phylo4-method nodeLabels
#' nodeLabels-methods nodeLabels,phylo4-method nodeLabels<-
#' nodeLabels<-,phylo4,character-method nodeLabels<-,phylo4d,ANY-method
#' tipLabels tipLabels-methods tipLabels,phylo4-method tipLabels<-
#' tipLabels<-,phylo4,character-method tipLabels<-,phylo4d,character-method
#' hasEdgeLabels hasEdgeLabels-methods hasEdgeLabels,phylo4-method edgeLabels
#' edgeLabels<- edgeLabels-methods edgeLabels,phylo4-method
#' edgeLabels<-,phylo4,character-method
#' @docType methods
#' @param x a phylo4 or phylo4d object.
#' @param object a phylo4 or phylo4d object.
#' @param type which type of labels: \code{all} (tips and internal nodes),
#' \code{tip} (tips only), \code{internal} (internal nodes only).
#' @param value a vector of class \code{character}, see Details for more
#' information.
#' @param use.names should the names of the vector used to create/update labels
#' be used to match the labels? See Details for more information.
#' @section Methods: \describe{ \item{labels}{\code{signature(object =
#' "phylo4")}: tip and/or internal node labels, ordered by node ID}
#' 
#' \item{hasDuplicatedLabels}{\code{signature(object = "phylo4")}: are any
#' labels duplicated?}
#' 
#' \item{tipLabels}{\code{signature(object = "phylo4")}: tip labels, ordered by
#' node ID}
#' 
#' \item{hasNodeLabels}{\code{signature(object = "phylo4")}: whether tree has
#' (internal) node labels} \item{nodeLabels}{\code{signature(object =
#' "phylo4")}: internal node labels, ordered by node ID}
#' 
#' \item{hasEdgeLabels}{\code{signature(object = "phylo4")}: whether tree has
#' (internal) edge labels} \item{edgeLabels}{\code{signature(object =
#' "phylo4")}: internal edge labels, ordered according to the edge matrix} }
#' @examples
#' 
#' 
#' data(geospiza)
#' 
#' ## Return labels from geospiza
#' tipLabels(geospiza)
#' 
#' ## Internal node labels in geospiza are empty
#' nodeLabels(geospiza)
#' 
#' ## Creating internal node labels
#' ndLbl <- paste("n", 1:nNodes(geospiza), sep="")
#' nodeLabels(geospiza) <- ndLbl
#' nodeLabels(geospiza)
#' 
#' ## naming the labels
#' names(ndLbl) <- nodeId(geospiza, "internal")
#' 
#' ## shuffling the labels
#' (ndLbl <- sample(ndLbl))
#' 
#' ## by default, the labels are attributed in the order
#' ## they are given:
#' nodeLabels(geospiza) <- ndLbl
#' nodeLabels(geospiza)
#' 
#' ## but use.names puts them in the correct order
#' labels(geospiza, "internal", use.names=TRUE) <- ndLbl
#' nodeLabels(geospiza)
#' 
## return labels in increasing node order
setMethod("labels", signature(object="phylo4"),
  function(object, type = c("all", "tip", "internal")) {
    type <- match.arg(type)
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    label <- object@label
    id <- nodeId(object, type)
    lbl <- label[match(id, names(label))]
    # reassign names b/c any unmatched will be NA (could instead assign
    # names only to the unmatched ones, but this seems simpler)
    names(lbl) <- id
    return(lbl)
})

setReplaceMethod("labels",
                 signature(x="phylo4", type="ANY",
                           use.names="ANY", value="character"),
   function(x, type = c("all", "tip", "internal"),
            use.names, ..., value) {

       ## Default options
       if(missing(type))
           type <- "all"
       if (missing(use.names))
           use.names <- FALSE

       type <- match.arg(type)

       ## generate new labels of the desired type
       new.label <- .createLabels(value, nTips(x), nNodes(x), use.names,
           type=type)

       ## replace existing labels and add new ones as needed
       old.label <- x@label
       old.index <- match(names(new.label), names(old.label))
       isNew <- is.na(old.index)
       old.label[old.index[!isNew]] <- new.label[!isNew]
       updated.label <- c(old.label, new.label[isNew])

       ## for efficiency, drop any NA labels
       x@label <- updated.label[!is.na(updated.label)]

       if(is.character(checkval <- checkPhylo4(x)))
           stop(checkval)
       else
           return(x)
   })

### Duplicated Labels
setMethod("hasDuplicatedLabels", signature(x="phylo4", type="ANY"),
  function(x, type=c("all", "tip", "internal")) {
      ## Default options
      if (missing(type)) {
          type <- "all"
      }
      type <- match.arg(type)
      any(duplicated(na.omit(labels(x, type))))
})


### Node Labels
setMethod("hasNodeLabels", signature(x="phylo4"),
 function(x) {
    !all(is.na(nodeLabels(x)))
})

setMethod("nodeLabels", signature(x="phylo4"),
 function(x) {
    labels(x, type="internal")
})

setReplaceMethod("nodeLabels", signature(x="phylo4", value="character"),
  function(x, ..., value) {
      labels(x, type="internal", ...) <- value
      if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
      x
  })

### Tip labels
setMethod("tipLabels", signature(x="phylo4"),
 function(x) {
    labels(x, type="tip")
    })

setReplaceMethod("tipLabels", signature(x="phylo4", value="character"),
  function(x, ...,  value) {
      labels(x, type="tip", ...) <- value
      if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
      return(x)
  })


### Edge labels
setMethod("hasEdgeLabels", signature(x="phylo4"),
 function(x) {
    !all(is.na(x@edge.label))
})

# return edge labels in order by edgeIds (same order as edge matrix)
setMethod("edgeLabels", signature(x="phylo4"),
  function(x) {
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    id <- edgeId(x, "all")
    lbl <- x@edge.label[match(id, names(x@edge.label))]
    names(lbl) <- id
    return(lbl)
})

setReplaceMethod("edgeLabels", signature(x="phylo4", value="character"),
  function(x, ..., value) {
    lbl <- .createEdge(value, x@edge, type="labels")    
    x@edge.label <- lbl[!is.na(lbl)]
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    x
  })


#########################################################
### Displaying phylo4: print, show, head, tail and summary
#########################################################

### print


#' print a phylogeny
#' 
#' Prints a phylo4 or phylo4d object in data.frame format with user-friendly
#' column names
#' 
#' This is a user-friendly version of the tree representation, useful for
#' checking that objects were read in completely and translated correctly. The
#' phylogenetic tree is represented as a list of numbered nodes, linked in a
#' particular way through time (or rates of evolutionary change).  The topology
#' is given by the pattern of links from each node to its ancestor. Also given
#' are the taxon names, node type (root/internal/tip) and phenotypic data (if
#' any) associated with the node, and the branch length from the node to its
#' ancestor. A list of nodes (descendants) and ancestors is minimally required
#' for a phylo4 object.
#' 
#' @param x a \code{phylo4} tree or \code{phylo4d} tree+data object
#' @param edgeOrder in the data frame returned, the option 'pretty' returns the
#' internal nodes followed by the tips, the option 'real' returns the nodes in
#' the order they are stored in the edge matrix.
#' @param printall default prints entire tree. printall=FALSE returns the first
#' 6 rows
#' @return A data.frame with a row for each node (descendant), sorted as
#' follows: root first, then other internal nodes, and finally tips.\cr The
#' returned data.frame has the following columns:\cr \item{label}{Label for the
#' taxon at the node (usually species name).} \item{node}{Node number, i.e. the
#' number identifying the node in \code{x@edge}.} \item{ancestor}{Node number
#' of the node's ancestor.} \item{branch.length}{The branch length connecting
#' the node to its ancestor (NAs if missing).} \item{node.type}{"root",
#' "internal", or "tip". (internally generated)} \item{data}{phenotypic data
#' associated with the nodes, with separate columns for each variable.}
#' @note This is the default show() method for phylo4, phylo4d. It prints the
#' user-supplied information for building a phylo4 object. For a full
#' description of the phylo4 S4 object and slots, see \code{\link{phylo4}}.
#' @author Marguerite Butler Thibaut Jombart
#' \email{jombart@@biomserv.univ-lyon1.fr} Steve Kembel
#' @keywords methods
#' @examples
#' 
#' 
#' tree.phylo <- ape::read.tree(text="((a,b),c);")
#' tree <- as(tree.phylo, "phylo4")
#' ##plot(tree,show.node=TRUE) ## plotting broken with empty node labels: FIXME
#' tip.data <- data.frame(size=c(1,2,3), row.names=c("a", "b", "c"))
#' treedata <- phylo4d(tree, tip.data)
#' plot(treedata)
#' print(treedata)
#' 
#' 
#' @export printphylo4
printphylo4 <- function(x, edgeOrder=c("pretty", "real"), printall=TRUE) {
    if(!nrow(edges(x))) {
        msg <- paste("Empty \'", class(x), "\' object\n", sep="")
        cat(msg)
    }
    else {
        toRet <- .phylo4ToDataFrame(x, edgeOrder)
        if (printall) {
            print(toRet)
        }
        else {
            print(head(toRet))
        }
    }
}

### Hack for print/show
### from http://tolstoy.newcastle.edu.au/R/e2/devel/06/12/1363.html
setMethod("print", "phylo4", printphylo4)
setMethod("show", signature(object="phylo4"),
   function(object) printphylo4(object))

### names
setMethod("names", signature(x="phylo4"),
 function(x) {
    temp <- rev(names(attributes(x)))[-1]
    return(rev(temp))
})

### Head and Tail
setMethod("head", signature(x="phylo4"),
  function(x, n=20) {
      head(as(x,"data.frame"),n=n)
  })

setMethod("tail", signature(x="phylo4"),
  function(x, n=20) {
      tail(as(x, "data.frame"), n=n)
  })

### summary
#' Displaying phylo4 object
#' 
#' Display methods for phylo4 and phylo4d phylogenetic trees
#' 
#' 
#' @name phylo4-display
#' @aliases print,phylo4-method show,phylo4-method head,phylo4-method
#' tail,phylo4-method summary,phylo4-method names,phylo4-method
#' @docType methods
#' @param x a phylo4 object
#' @param object a phylo4 object
#' @param edgeOrder Character string indicating whether the edges should be
#' printed as ordered in the tree "real" (e.g. preorder or postorder), or
#' "pretty" printed with tips collated together
#' @param printall If TRUE all tip labels are printed
#' @param quiet a logical stating whether the results of the summary should be
#' printed to the screen (FALSE, default) or not (TRUE)
#' @return
#' 
#' The \code{summary} method invisibly returns a list with the following
#' components:
#' 
#' \item{list("name")}{the name of the object} \item{list("nb.tips")}{the
#' number of tips} \item{list("nb.nodes")}{the number of nodes}
#' \item{list("mean.el")}{mean of edge lengths} \item{list("var.el")}{variance
#' of edge lengths (estimate for population)} \item{list("sumry.el")}{summary
#' (i.e. range and quartiles) of the edge lengths}
#' \item{list("degree")}{(optional) degree (i.e. number of descendants) of each
#' node; displayed only when there are polytomies}
#' \item{list("polytomy")}{(optional) type of polytomy for each node:
#' \sQuote{node}, \sQuote{terminal} (all descendants are tips) or
#' \sQuote{internal} (at least one descendant is an internal node); displayed
#' only when there are polytomies}
#' 
#' The \code{names} method returns a vector of characters corresponding to the
#' names of the slots.
#' @section Methods: \describe{ \item{print}{\code{signature(x = "phylo4")}:
#' print method} \item{show}{\code{signature(object = "phylo4")}: show method }
#' \item{summary}{\code{signature(object = "phylo4")}: summary method}
#' \item{names}{\code{signature(x = "phylo4")}: gives the slot names}
#' \item{head}{\code{signature(object = "phylo4")}: show first few nodes}
#' \item{tail}{\code{signature(object = "phylo4")}: show last few nodes} }
#' @author Ben Bolker, Thibaut Jombart
#' @seealso The \code{\link{phylo4}} constructor, the \code{\link{checkPhylo4}}
#' function to check the validity of \code{phylo4} objects. See also the
#' \code{\link{phylo4d}} constructor and the \linkS4class{phylo4d} class.
#' @keywords methods
#' @examples
#' 
#' 
#'   tOwls <- "(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
#'   tree.owls <- ape::read.tree(text=tOwls)
#'   P1 <- as(tree.owls, "phylo4")
#'   P1
#'   summary(P1)
#' 
#' 
#'   ## summary of a polytomous tree
#'   E <- matrix(c(
#'       8,  9,
#'       9, 10,
#'      10,  1,
#'      10,  2,
#'       9,  3,
#'       9,  4,
#'       8, 11,
#'      11,  5,
#'      11,  6,
#'      11,  7,
#'       0,  8), ncol=2, byrow=TRUE)
#' 
#'   P2 <- phylo4(E)
#'   nodeLabels(P2) <- as.character(nodeId(P2, "internal"))
#'   plot(P2, show.node.label=TRUE)
#'   sumryP2 <- summary(P2)
#'   sumryP2
#' 
setMethod("summary", signature(object="phylo4"),
  function(object, quiet=FALSE) {

      res <- list()

      ## build the result object
      res$name <- deparse(substitute(object, sys.frame(-1)))
      res$nb.tips <- nTips(object)
      res$nb.nodes <- nNodes(object)

      if(hasEdgeLength(object)) {
          edge.length <- edgeLength(object)
          res$mean.el <- mean(edge.length, na.rm=TRUE)
          res$var.el <- var(edge.length, na.rm=TRUE)
          if (isRooted(object) && is.na(edgeLength(object, rootNode(object)))) {
              root.index <- match(edgeId(object, "root"), names(edge.length))
              res$sumry.el <- summary(edge.length[-root.index])
          } else {
              res$sumry.el <- summary(edge.length)
          }
      }

      ## check for polytomies
      if (hasPoly(object)) {
          E <- edges(object)
          temp <- tabulate(na.omit(E[,1]))
          degree <- temp[na.omit(E[,1])] # contains the degree of the ancestor for all edges
          endsAtATip <- !(E[,2] %in% E[,1])
          terminPoly <- (degree>2) & endsAtATip
          internPoly <- (degree>2) & !endsAtATip
          res$degree <- degree
          res$polytomy <- rep("none",nrow(E))
          res$polytomy[terminPoly] <- "terminal"
          res$polytomy[internPoly] <- "internal"
          ## now just keep information about nodes (not all edges)
          nod <- unique(E[,1])
          idx <- match(nod,E[,1])
          res$degree <- res$degree[idx]
          names(res$degree) <- nodeLabels(object)
          res$polytomy <- res$polytomy[idx]
          names(res$polytomy) <- nodeLabels(object)
      }

      ## model info
      res$loglik <- attr(object, "loglik")
      res$para <- attr(object, "para")
      res$xi <- attr(object, "xi")

      ## if quiet, stop here
      if(quiet) return(invisible(res))

      ## now, print to screen is !quiet
      cat("\n Phylogenetic tree :", res$name, "\n\n")
      cat(" Number of tips    :", res$nb.tips, "\n")
      cat(" Number of nodes   :", res$nb.nodes, "\n")
      ## cat("  ")
      if(hasEdgeLength(object)) {
          cat(" Branch lengths:\n")
          cat("        mean         :", res$mean.el, "\n")
          cat("        variance     :", res$var.el, "\n")
          cat("        distribution :\n")
          print(res$sumry.el)
      }
      else {
          cat(" Branch lengths    : No branch lengths.\n")
      }
      if (hasPoly(object)) {
          cat("\nDegree of the nodes  :\n")
          print(res$degree)
          cat("\n")
          cat("Types of polytomy:\n")
          print(res$polytomy)
          cat("\n")
      }

      if (!is.null(attr(object, "loglik"))) {
          cat("Phylogeny estimated by maximum likelihood.\n")
          cat("  log-likelihood:", attr(object, "loglik"), "\n\n")
          npart <- length(attr(object, "para"))
          for (i in 1:npart) {
              cat("partition ", i, ":\n", sep = "")
              print(attr(object, "para")[[i]])
              if (i == 1)
                  next
              else cat("  contrast parameter (xi):", attr(object,"xi")[i - 1], "\n")
        }
      }
      return(invisible(res))

  }) # end setMethod summary phylo4




#########################################################
### Ordering
#########################################################

#' reordering trees within phylobase objects
#' 
#' Methods for reordering trees into various traversal orders
#' 
#' The \code{reorder} method takes a \code{phylo4} or \code{phylo4d} tree and
#' orders the edge matrix (i.e. \code{edges(x)}) in the requested traversal
#' order. Currently only two orderings are permitted, and both require rooted
#' trees. In "postorder", a node's descendants come before that node, thus the
#' root, which is ancestral to all nodes, comes last.  In "preorder", a node is
#' visited before its descendants, thus the root comes first.
#' 
#' A method is also defined that takes an \code{ape phylo} object.  This also
#' takes an order argument, however, 'pruningwise' and 'cladewise' are the only
#' acceptable parameters. This is because this method actually uses the
#' \code{ape reorder()} command to complete the ordering.
#' 
#' @name reorder-methods
#' @aliases reorder-methods reorder,phylo-method reorder,phylo4-method
#' reorder,phylo4d-method
#' @docType methods
#' @param x a \code{phylo4} or \code{phylo4d} object
#' @param order The desired traversal order; currently only 'preorder' and
#' 'postorder' are allowed for \code{phylo4} and \code{phylo4d} objects,
#' whereas only 'cladewise' and 'pruningwise' are allowed for \code{phylo}
#' objects
#' @return A \code{phylo4} or \code{phylo4d} object with the edge, label,
#' length and data slots ordered as \code{order}, which is itself recorded in
#' the order slot.
#' @note The "preorder" parameter corresponds to "cladewise" in the \code{ape}
#' package, and "postorder" corresponds (almost but close enough?) to
#' "pruningwise".
#' 
#' See \url{http://ape.mpl.ird.fr/misc/FormatTreeR_28July2008.pdf}
#' @section Methods: \describe{ \item{x = "phylo"}{reorders a \code{phylo}
#' object} \item{x = "phylo4"}{reorders a \linkS4class{phylo4} object} \item{x
#' = "phylo4d"}{reorders a \linkS4class{phylo4d} object} }
#' @author Peter Cowan, Jim Regetz
#' @seealso \code{\link[ape]{reorder.phylo}} in the \code{ape} package.
#' \code{\link{ancestors}} \code{\link{ancestor}} \code{\link{siblings}}
#' \code{\link{children}} \code{\link{descendants}}
#' @keywords methods
#' @examples
#' phy <- phylo4(ape::rtree(5))
#' edges(reorder(phy, "preorder"))
#' edges(reorder(phy, "postorder"))
orderIndex <- function(x, order=c("preorder", "postorder")) {

    order <- match.arg(order)
    if(!isRooted(x)){
        stop("Tree must be rooted to reorder")
    }
    ## get a root node free edge matrix
    edge <- edges(x, drop.root=TRUE)
    ## Sort edges -- ensures that starting order of edge matrix doesn't
    ## affect the order of reordered trees
    edge <- edge[order(edge[, 2]), ]

    # recast order argument as integer to pass to C
    if(order == 'postorder') {
        iOrder <- 0L
    } else if(order == 'preorder') {
        iOrder <- 1L
    } else {stop(paste("Method for", order, "not implemented"))}

    if (!hasPoly(x) & !hasSingle(x)) {
        # method 1: faster, but only works if all internal nodes have
        # exactly two children (true binary tree)

        # extract nodes, separating descendants into left (first
        # encountered) and right (second encountered) for each ancestor
        isFirst <- !duplicated(edge[, 1])
        ancestor <- as.integer(edge[isFirst, 1])
        left <- as.integer(edge[isFirst, 2])
        right <- as.integer(edge[!isFirst, 2])[match(ancestor,
            edge[!isFirst, 1])]
        descendantNew <- rep(0L, nEdges(x))
        root <- as.integer(rootNode(x))
        nEdge <- as.integer(length(ancestor))

        descendantReord <- .C("reorderBinary", descendantNew, root,
            ancestor, left, right, nEdge, iOrder)[[1]]

    } else {
        # method 2: not as fast, but robust to singletons and polytomies

        # extract ancestors and descendants
        ancestor <- as.integer(edge[,1])
        descendant <- as.integer(edge[,2])
        descendantNew <- rep(0L, nEdges(x))
        root <- as.integer(rootNode(x))
        nEdge <- as.integer(nrow(edge))

        descendantReord <- .C("reorderRobust", descendantNew, root,
            ancestor, descendant, nEdge, iOrder)[[1]]

    }

    ## Original pure R implementation of the above:
    #### recursive functions are placed first and calls to those functions below
    ##postOrder <- function(node) {
    ##    ## this function returns a vector of nodes in the post order traversal
    ##    ## get the descendants
    ##    traversal <- NULL
    ##    ## edge -- defined above, outside this function
    ##    ## extensive testing found this loop to be faster than apply() etc.
    ##    for(i in edge[edge[, 1] == node, 2]) {
    ##        traversal <- c(traversal, postOrder(i))
    ##    }
    ##    c(traversal, node)
    ##}
    ##preOrder  <- function(node) {
    ##    ## see expanded code in comments of postOrder()
    ##    ## only difference here is that we record current node, then descendants
    ##    traversal <- NULL
    ##    for(i in edge[edge[, 1] == node, 2]) {
    ##        traversal <- c(traversal, preOrder(i))
    ##    }
    ##    c(node, traversal)
    ##}
    ##if(order == 'postorder') {
    ##    descendantReord <- postOrder(rootNode(x))
    ##} else if(order == 'preorder') {
    ##    descendantReord <- preOrder(rootNode(x))
    ##} else {stop(paste("Method for", order, "not implemented"))}

    ## match the new node order to the old order to get an index
    index <- match(descendantReord, edges(x)[, 2])

}

setMethod("reorder", signature(x="phylo4"),
 function(x, order=c("preorder", "postorder")) {
    ## call orderIndex and use that index to order edges, labels and lengths
    order   <- match.arg(order)
    index   <- orderIndex(x, order)
    x@order <- order
    x@edge  <- edges(x)[index, ]
    if(hasEdgeLabels(x)) {
        x@edge.label  <- x@edge.label[index]
    }
    if(hasEdgeLength(x)) {
        x@edge.length <- x@edge.length[index]
    }
    x
})


#########################################################
### Tree properties
#########################################################

setMethod("isUltrametric", signature(x="phylo4"),
  function(x, tol=.Machine$double.eps^.5) {
    if (!hasEdgeLength(x)) {
      stop("The tree has no edge lengths.")
    }
    if (identical(all.equal.numeric(var(depthTips(x)), 0, tolerance=tol), TRUE)) {
      TRUE
    }
    else FALSE
  })
