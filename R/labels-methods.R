
#########################################################
### Label accessors
#########################################################

# @aliases labels<- labels,phylo4-method
# labels<-,phylo4,ANY,ANY,character-method
# labels<-,phylo4d,ANY,ANY,character-method hasDuplicatedLabels
# hasDuplicatedLabels-methods hasDuplicatedLabels,phylo4-method hasNodeLabels
# hasNodeLabels-methods hasNodeLabels,phylo4-method nodeLabels
# nodeLabels-methods nodeLabels,phylo4-method nodeLabels<-
# nodeLabels<-,phylo4,character-method nodeLabels<-,phylo4d,ANY-method
# tipLabels tipLabels-methods tipLabels,phylo4-method tipLabels<-
# tipLabels<-,phylo4,character-method tipLabels<-,phylo4d,character-method
# hasEdgeLabels hasEdgeLabels-methods hasEdgeLabels,phylo4-method edgeLabels
# edgeLabels<- edgeLabels-methods edgeLabels,phylo4-method
# edgeLabels<-,phylo4,character-method

##' Labels for phylo4/phylo4d objects
##' 
##' Methods for creating, accessing and updating labels in
##' phylo4/phylo4d objects
##' 
##' In phylo4/phylo4d objects, tips must have labels (that's why there
##' is no method for hasTipLabels), internal nodes and edges can have
##' labels.
##' 
##' Labels must be provided as a vector of class \code{character}. The
##' length of the vector must match the number of elements they label.
##' 
##' The option \code{use.names} allows the user to match a label to a
##' particular node. In this case, the vector must have names that
##' match the node numbers.
##' 
##' The function \code{labels} is mostly intended to be used
##' internally.
##' 
##' @name phylo4-labels
##' @aliases labels
##' @docType methods
##' @param x a phylo4 or phylo4d object.
##' @param object a phylo4 or phylo4d object.
##' @param type which type of labels: \code{all} (tips and internal nodes),
##' \code{tip} (tips only), \code{internal} (internal nodes only).
##' @param value a vector of class \code{character}, see Details for more
##' information.
##' @param use.names should the names of the vector used to create/update labels
##' be used to match the labels? See Details for more information.
##' @section Methods: \describe{ \item{labels}{\code{signature(object =
##' "phylo4")}: tip and/or internal node labels, ordered by node ID}
##' 
##' \item{hasDuplicatedLabels}{\code{signature(object = "phylo4")}: are any
##' labels duplicated?}
##' 
##' \item{tipLabels}{\code{signature(object = "phylo4")}: tip labels, ordered by
##' node ID}
##' 
##' \item{hasNodeLabels}{\code{signature(object = "phylo4")}: whether tree has
##' (internal) node labels} \item{nodeLabels}{\code{signature(object =
##' "phylo4")}: internal node labels, ordered by node ID}
##' 
##' \item{hasEdgeLabels}{\code{signature(object = "phylo4")}: whether tree has
##' (internal) edge labels} \item{edgeLabels}{\code{signature(object =
##' "phylo4")}: internal edge labels, ordered according to the edge matrix} }
##' @export
##' @rdname labels-methods
##' @include phylo4-class.R phylo4-methods.R phylo4-accessors.R nodeId-methods.R
##' @author Ben Bolker, Peter Cowan, Steve Kembel, Francois Michonneau
##' @return labels in ascending order.
##' @examples
##' 
##' data(geospiza)
##' 
##' ## Return labels from geospiza
##' tipLabels(geospiza)
##' 
##' ## Internal node labels in geospiza are empty
##' nodeLabels(geospiza)
##' 
##' ## Creating internal node labels
##' ndLbl <- paste("n", 1:nNodes(geospiza), sep="")
##' nodeLabels(geospiza) <- ndLbl
##' nodeLabels(geospiza)
##' 
##' ## naming the labels
##' names(ndLbl) <- nodeId(geospiza, "internal")
##' 
##' ## shuffling the labels
##' (ndLbl <- sample(ndLbl))
##' 
##' ## by default, the labels are attributed in the order
##' ## they are given:
##' nodeLabels(geospiza) <- ndLbl
##' nodeLabels(geospiza)
##' 
##' ## but use.names puts them in the correct order
##' labels(geospiza, "internal", use.names=TRUE) <- ndLbl
##' nodeLabels(geospiza)
setGeneric("labels")

##' @rdname labels-methods
##' @aliases labels,phylo4-method
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

##' @rdname labels-methods
##' @aliases labels<-
setGeneric("labels<-",
           function(x, type, use.names, ..., value) {
               standardGeneric("labels<-")
           })

##' @name labels<-
##' @rdname labels-methods
##' @aliases labels<-,phylo4,ANY,ANY,character-method
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

##### -------- hasDuplicatedLabels

##' @rdname labels-methods
##' @aliases hasDuplicatedLabels
setGeneric("hasDuplicatedLabels",
           function(x, type) {
               standardGeneric("hasDuplicatedLabels")
           })

##' @rdname labels-methods
##' @aliases hasDuplicatedLabels,phylo4,ANY-method
setMethod("hasDuplicatedLabels", signature(x="phylo4", type="ANY"),
  function(x, type=c("all", "tip", "internal")) {
      ## Default options
      if (missing(type)) {
          type <- "all"
      }
      type <- match.arg(type)
      any(duplicated(na.omit(labels(x, type))))
})

##### --------- hasNodeLabels

##' @rdname labels-methods
##' @aliases hasNodeLabels
setGeneric("hasNodeLabels", function(x) {
    standardGeneric("hasNodeLabels")
})

##' @rdname labels-methods
##' @aliases hasNodeLabels,phylo4-method
setMethod("hasNodeLabels", signature(x="phylo4"),
 function(x) {
    !all(is.na(nodeLabels(x)))
})

##### ---------- nodeLabels

##' @rdname labels-methods
##' @aliases nodeLabels
setGeneric("nodeLabels", function(x) {
    standardGeneric("nodeLabels")
})

##' @rdname labels-methods
##' @aliases nodeLabels,phylo4-method
setMethod("nodeLabels", signature(x="phylo4"),
 function(x) {
    labels(x, type="internal")
})

##' @rdname labels-methods
##' @aliases nodeLabels<-
setGeneric("nodeLabels<-",
           function(x, ..., value) {
               standardGeneric("nodeLabels<-")
           })

##' @name nodeLabels<-
##' @rdname labels-methods
##' @aliases nodeLabels<-,phylo4,character-method
setReplaceMethod("nodeLabels", signature(x="phylo4", value="character"),
  function(x, ..., value) {
      labels(x, type="internal", ...) <- value
      if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
      x
  })

##### ---------- tipLabels

##' @rdname labels-methods
##' @aliases tipLabels
setGeneric("tipLabels", function(x) {
    standardGeneric("tipLabels")
})

##' @rdname labels-methods
##' @aliases tipLabels,phylo4-method
setMethod("tipLabels", signature(x="phylo4"),
 function(x) {
    labels(x, type="tip")
    })

##' @rdname labels-methods
##' @aliases tipLabels<-
setGeneric("tipLabels<-",
   function(x, ..., value) {
       standardGeneric("tipLabels<-")
})

##' @name tipLabels<-
##' @rdname labels-methods
##' @aliases tipLabels<-,phylo4,character-method
setReplaceMethod("tipLabels", signature(x="phylo4", value="character"),
  function(x, ...,  value) {
      labels(x, type="tip", ...) <- value
      if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
     return(x)
  })


##### ---------- hasEdgeLabels

##' @rdname labels-methods
##' @aliases hasEdgeLabels
setGeneric("hasEdgeLabels", function(x) {
    standardGeneric("hasEdgeLabels")
})

##' @rdname labels-methods
##' @aliases hasEdgeLabels,phylo4-method
setMethod("hasEdgeLabels", signature(x="phylo4"),
 function(x) {
    !all(is.na(x@edge.label))
})

##### ----------  edgeLabels

##' @rdname labels-methods
##' @aliases edgeLabels
setGeneric("edgeLabels", function(x) {
    standardGeneric("edgeLabels")
})

##' @rdname labels-methods
##' @aliases edgeLabels,phylo4-method
setMethod("edgeLabels", signature(x="phylo4"),
  function(x) {
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    id <- edgeId(x, "all")
    lbl <- x@edge.label[match(id, names(x@edge.label))]
    names(lbl) <- id
    return(lbl)
})

##' @rdname labels-methods
##' @aliases edgeLabels<-
setGeneric("edgeLabels<-",
           function(x, ..., value) {
               standardGeneric("edgeLabels<-")
           })

##' @name edgeLabels<-
##' @rdname labels-methods
##' @aliases edgeLabels<-,phylo4,character-method
setReplaceMethod("edgeLabels", signature(x="phylo4", value="character"),
  function(x, ..., value) {
    lbl <- .createEdge(value, x@edge, type="labels")    
    x@edge.label <- lbl[!is.na(lbl)]
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    x
  })

