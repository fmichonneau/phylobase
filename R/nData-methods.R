
##' Retrieves the number of datasets in phylo4d objects
##' 
##' Method to retrieve the number of datasets associated with a phylogenetic
##' tree stored as a phylo4d object
##' 
##' \code{nData} returns the number of datasets (i.e., columns) that are
##' associated with a \code{phylo4d} object.
##' 
##' @param x A \code{phylo4d} object
##' @return \code{nData} returns a vector.
##' @author Francois Michonnea
##' @seealso \code{\link{tdata}}, \code{\link{phylo4d}}
##' @keywords methods
##' @export
##' @rdname nData-methods
##' @docType methods
##' @examples
##' 
##'   data(geospiza)
##'   nData(geospiza)
setGeneric("nData", function(x, ...) {
     standardGeneric("nData")
})

##' @rdname nData-methods
##' @aliases nData,phylo4d-method
setMethod("nData", signature(x="phylo4d"), function(x) {
    ncol(x@data)
})
