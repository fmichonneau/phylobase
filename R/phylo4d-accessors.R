
##' Tests for presence of data associated with trees stored as phylo4d objects
##' 
##' Methods that test for the presence of data associated with trees stored as
##' phylo4d objects.
##' 
##' The outcome of the test is based on row names of the data frame stored in
##' \code{data}. If there are no rows having row names from the set
##' \code{nodeId(x, "tip")}, then \code{hasTipData} returns FALSE.  Likewise, if
##' there are no rows having row names from the set \code{nodeId(x,
##' "internal")}, then \code{hasNodeData} returns FALSE.
##' 
##' @aliases hasNodeData hasNodeData-methods hasNodeData,phylo4d-method
##' hasTipData hasTipData-methods hasTipData,phylo4d-method
##' @param x a phylo4d object
##' @return \item{list("logical")}{return \code{TRUE} or \code{FALSE} depending
##' whether data are associated with the tree (i.e., the slots \code{tip.data}
##' or \code{node.data} are not empty)}
##' @section Methods: \describe{ \item{hasNodeData}{\code{signature(object =
##' "phylo4d")}: whether tree has internal node data}
##' \item{hasTipData}{\code{signature(object = "phylo4d")}: whether tree has
##' data associated with its tips} }
##' @author Ben Bolker, Thibault Jombart, Francois Michonneau
##' @seealso \code{\link{phylo4d}} constructor and \code{\linkS4class{phylo4d}}
##' class.
##' @name hasTipData
##' @keywords methods
##' @docType methods
##' @include phylo4d-class.R
##' @export
##' @examples
##'   data(geospiza)
##'   hasTipData(geospiza)  ## TRUE
##'   hasNodeData(geospiza) ## FALSE
##'
setGeneric("hasTipData", function(x) {
    standardGeneric("hasTipData")
})

##' @name hasTipData-methods
##' @rdname hasTipData
##' @aliases hasTipData-method,phylo4d-method
setMethod("hasTipData", signature(x="phylo4d"),
 function(x) {
    ncol(tdata(x, type="tip", empty.columns=FALSE)) > 0
})

##' @rdname hasTipData
##' @aliases hasNodeData-methods
setGeneric("hasNodeData", function(x) {
    standardGeneric("hasNodeData")
})

##' @name hasNodeData-methods
##' @rdname hasTipData
##' @aliases hasNodeData,phylo4d-method
setMethod("hasNodeData", signature(x="phylo4d"),
 function(x) {
    ncol(tdata(x, type="internal", empty.columns=FALSE)) > 0
})
