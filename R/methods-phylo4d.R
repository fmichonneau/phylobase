
#' Retrieving or updating tip and node data in phylo4d objects
#' 
#' Methods to retrieve or update tip, node or all data associated with a
#' phylogenetic tree stored as a phylo4d object
#' 
#' 
#' @aliases tdata tdata-method tdata,phylo4d-method tdata<-
#' tdata<-,phylo4d-method tdata<-,phylo4d,ANY-method tipData tipData-method
#' tipData,phylo4d-method tipData<- tipData<-,phylo4d-method
#' tipData<-,phylo4d,ANY-method nodeData nodeData-method
#' nodeData,phylo4d-method nodeData<- nodeData<-,phylo4d-method
#' nodeData<-,phylo4d,ANY-method
#' @param x A \code{phylo4d} object
#' @param type The type of data to retrieve or update: \dQuote{\code{all}}
#' (default) for data associated with both tip and internal nodes,
#' \dQuote{\code{tip}} for data associated with tips only,
#' \dQuote{\code{internal}} for data associated with internal nodes only.
#' @param label.type How should the tip/node labels from the tree be returned?
#' \dQuote{\code{row.names}} returns them as row names of the data frame,
#' \dQuote{\code{column}} returns them in the first column of the data frame.
#' This options is useful in the case of missing (\code{NA}) or non-unique
#' labels.
#' @param empty.columns Should columns filled with \code{NA} be returned?
#' @param merge.data if tip or internal node data are provided and data already
#' exists for the other type, this determines whether columns with common names
#' will be merged together (default TRUE). If FALSE, columns with common names
#' will be preserved separately, with \dQuote{.tip} and \dQuote{.node} appended
#' to the names. This argument has no effect if tip and node data have no
#' column names in common, or if type=\dQuote{all}.
#' @param clear.all If only tip or internal node data are to be replaced,
#' should data of the other type be dropped?
#' @param \dots For the tipData and nodeData accessors, further arguments to be
#' used by tdata. For the replacement forms, further arguments to be used by
#' \code{formatData} (e.g.  \code{match.data}), see \link{formatData} for more
#' details.
#' @param value a data frame (or object to be coerced to one) to replace the
#' values associated with the nodes specified by the argument \code{type}
#' @return \code{tdata} returns a data frame
#' @section Methods: \describe{
#' \item{tdata}{\code{signature(object="phylo4d")}: retrieve or update data
#' associated with a tree in a \code{phylo4d} object} }
#' @author Ben Bolker, Thibaut Jombart, Francois Michonneau
#' @seealso \code{\link{phylo4d}}
#' @keywords methods
#' @examples
#' 
#'    data(geospiza)
#'    tdata(geospiza)
#'    tipData(geospiza) <- 1:nTips(geospiza)
#'    tdata(geospiza)
#'    \dontshow{data(geospiza)}
#' 
setMethod("tdata", signature(x="phylo4d"),
  function(x, type=c("all", "tip", "internal"),
           label.type=c("row.names","column"),
           empty.columns=TRUE) {

      ## Returns data associated with the tree
      ## Note: the function checks for unique labels. It's currently unecessary
      ## but could be useful in the future if non-unique labels are allowed.

      type <- match.arg(type)
      label.type <- match.arg(label.type)

      ids <- nodeId(x, type)
      labs <- labels(x, type)
      ## replace any missing labels with node numbers
      labs[is.na(labs)] <- names(labs)[is.na(labs)]
      tdata <- x@data[match(ids, row.names(x@data)), , drop=FALSE]
      row.names(tdata) <- ids
      data.names <- labs[match(names(labs), rownames(tdata))]

      if (label.type == "row.names") {
          if (!any(duplicated(data.names)) &&
              ## length(data.names) > 0 &&
              !any(is.na(data.names)) ) {
              row.names(tdata) <- data.names
          }
          else {
              warning("Non-unique or missing labels found, ",
                      "labels cannot be coerced to tdata row.names. ",
                      "Use the label.type argument to include labels ",
                      "as first column of data.")
          }
      }
      if (identical(label.type,"column")) {
          tdata <- data.frame(label=data.names, tdata)
      }

      ## remove empty columns (filled with NAs)
      if(!empty.columns) {
          emptyCol <- apply(tdata, 2, function(x) all(is.na(x)))
          tdata <- tdata[, !emptyCol, drop=FALSE]
      }

      tdata
  })

setReplaceMethod("tdata", signature(x="phylo4d", value="ANY"),
    function(x, type = c("all", "tip", "internal"), merge.data = TRUE,
        clear.all = FALSE, ..., value) {

    type <- match.arg(type)

    ## format new data
    value <- formatData(x, value, type, keep.all=TRUE, ...)

    ## get old data to keep (if any)
    if (clear.all || type=="all") {
        keep <- NULL
    } else {
        if (type=="tip") {
            keep <- tdata(x, type="internal", empty.column=FALSE)
            keep <- formatData(x, keep, "internal", match.data=FALSE)
        } else if (type=="internal") {
            keep <- tdata(x, type="tip", empty.column=FALSE)
            keep <- formatData(x, keep, "tip", match.data=FALSE)
        }
    }

    ## create updated data
    updated.data <- switch(type,
        tip = .phylo4Data(x, tip.data=value, node.data=keep,
            merge.data=merge.data),
        internal = .phylo4Data(x, tip.data=keep, node.data=value,
            merge.data=merge.data),
        all = .phylo4Data(x, all.data=value, merge.data=merge.data))

    ## try to arrange new columns after old columns
    kept <- names(updated.data) %in% names(keep)
    old.cols <- names(updated.data)[kept]
    new.cols <- names(updated.data)[!kept]
    x@data <- updated.data[c(old.cols, new.cols)]

    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    return(x)
})

### Tip data wrappers
setMethod("tipData", signature(x="phylo4d"), function(x, ...) {
    tdata(x, type="tip", ...)
})

setReplaceMethod("tipData", signature(x="phylo4d", value="ANY"),
    function(x, ...,  value) {
    tdata(x, type="tip", ...) <- value
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    return(x)
})

### Node data wrappers
setMethod("nodeData", signature(x="phylo4d"), function(x, ...) {
    tdata(x, type="internal", ...)
})

setReplaceMethod("nodeData", signature(x="phylo4d", value="ANY"),
    function(x, ...,  value) {
    tdata(x, type="internal", ...) <- value
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    return(x)
})

#' Adding data to a phylo4 or a phylo4d object
#' 
#' \code{addData} adds data to a \code{phylo4} (converting it in a
#' \code{phylo4d} object) or to a \code{phylo4d} object
#' 
#' Rules for matching data to tree nodes are identical to those used by the
#' \code{\link{phylo4d}} constructor.
#' 
#' If any column names in the original data are the same as columns in the new
#' data, ".old" is appended to the former column names and ".new" is appended
#' to the new column names.
#' 
#' The option \code{pos} is ignored (silently) if \code{x} is a \code{phylo4}
#' object. It is provided for compatibility reasons.
#' 
#' @aliases addData addData-methods addData,phylo4-method
#' addData,phylo4d-method
#' @param x a phylo4 or a phylo4d object
#' @param tip.data a data frame (or object to be coerced to one) containing
#' only tip data
#' @param node.data a data frame (or object to be coerced to one) containing
#' only node data
#' @param all.data a data frame (or object to be coerced to one) containing
#' both tip and node data
#' @param merge.data if both \code{tip.data} and \code{node.data} are provided,
#' it determines whether columns with common names will be merged together
#' (default TRUE). If FALSE, columns with common names will be preserved
#' separately, with ".tip" and ".node" appended to the names. This argument has
#' no effect if \code{tip.data} and \code{node.data} have no column names in
#' common.
#' @param pos should the new data provided be bound \code{before} or
#' \code{after} the pre-existing data?
#' @param \dots additional arguments to be passed to \link{formatData}
#' @return \code{addData} returns a \code{phylo4d} object.
#' @author Francois Michonneau
#' @seealso \code{\link{tdata}} for extracting or updating data and
#' \code{\link{phylo4d}} constructor.
#' @keywords methods
#' @examples
#' 
#'   data(geospiza)
#'   nDt <- data.frame(a=rnorm(nNodes(geospiza)), b=1:nNodes(geospiza),
#'     row.names=nodeId(geospiza, "internal"))
#'   t1 <- addData(geospiza, node.data=nDt)
setMethod("addData", signature(x="phylo4d"),
  function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
           merge.data=TRUE, pos=c("after", "before"), ...) {

    pos <- match.arg(pos)

    ## apply formatData to ensure data have node number rownames and
    ## correct dimensions
    tip.data <- formatData(phy=x, dt=tip.data, type="tip", ...)
    node.data <- formatData(phy=x, dt=node.data, type="internal", ...)
    all.data <- formatData(phy=x, dt=all.data, type="all", ...)
    ## combine data as needed
    new.data <- .phylo4Data(x=x, tip.data=tip.data, node.data=node.data,
        all.data=all.data, merge.data=merge.data)

    if (all(dim(new.data) == 0)) {
        return(x)
    }
    if (all(dim(x@data) == 0)) {
        x@data <- new.data
        return(x)
    }

    if (identical(pos, "after")) {
        new.data <- merge(x@data, new.data, by=0, all=TRUE,
            sort=FALSE, suffixes=c(".old", ".new"))
    } else {
        new.data <- merge(new.data, x@data, by=0, all=TRUE,
            sort=FALSE, suffixes=c(".new", ".old"))
    }
    row.names(new.data) <- new.data[["Row.names"]]
    x@data <- subset(new.data, select=-Row.names)

    x
})

setMethod("addData", signature(x="phylo4"),
  function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
           merge.data=TRUE, pos=c("after", "before"), ...) {
    phylo4d(x, tip.data=tip.data, node.data=node.data, all.data=all.data,
            merge.data=merge.data, ...)
})


#' Retrieves the number of datasets in phylo4d objects
#' 
#' Method to retrieve the number of datasets associated with a phylogenetic
#' tree stored as a phylo4d object
#' 
#' \code{nData} returns the number of datasets (i.e., columns) that are
#' associated with a \code{phylo4d} object.
#' 
#' @aliases nData nData,phylo4d-method
#' @param x A \code{phylo4d} object
#' @return \code{nData} returns a vector.
#' @author Francois Michonnea
#' @seealso \code{\link{tdata}}, \code{\link{phylo4d}}
#' @keywords methods
#' @examples
#' 
#'   data(geospiza)
#'   nData(geospiza)
setMethod("nData", signature(x="phylo4d"), function(x) {
    ncol(x@data)
})


#' Display methods for phylo4d objects
#' 
#' Methods used to display information about the data and the tree for phylo4d
#' objects.
#' 
#' 
#' @name phylo4d-display
#' @aliases show,phylo4d-method print,phylo4d-method head,phylo4d-method
#' tail,phylo4d-method summary,phylo4d-method names,phylo4d-method
#' @docType methods
#' @param x a phylo4d object
#' @param object a phylo4d object
#' @param edgeOrder Character string indicating whether the edges should be
#' printed as ordered in the tree "real" (e.g. preorder or postorder), or
#' "pretty" printed with tips collated together
#' @param quiet Should the summary be displayed on screen?
#' @param printall If TRUE all tip labels are printed
#' @return The \code{summary} method invisibly returns a list with the
#' following components: \item{list("name")}{the name of the object}
#' \item{list("nb.tips")}{the number of tips} \item{list("nb.nodes")}{the
#' number of nodes} \item{list("mean.el")}{mean of edge lengths}
#' \item{list("var.el")}{variance of edge lengths (estimate for population) }
#' \item{list("sumry.el")}{summary (i.e. range and quartiles) of the edge
#' lengths} \item{list("degree")}{(optional) type of polytomy for each node:
#' \sQuote{node}, \sQuote{terminal} (all descendants are tips) or
#' \sQuote{internal} (at least one descendant is an internal node); displayed
#' only when there are polytomies} \item{list("sumry.tips")}{(optional) summary
#' for the data associated with the tips} \item{list("sumry.nodes")}{(optional)
#' summary for the data associated with the internal nodes}
#' 
#' The \code{names} method returns a vector of characters corresponding to the
#' names of the slots.
#' @section Methods: \describe{ \item{print}{\code{signature(x = "phylo4d")}:
#' print method} \item{show}{\code{signature(object = "phylo4d")}: show method
#' } \item{summary}{\code{signature(object = "phylo4d")}: summary method}
#' \item{names}{\code{signature(x = "phylo4d")}: gives the slots names}
#' \item{head}{\code{signature(object = "phylo4d")}: show first few nodes}
#' \item{tail}{\code{signature(object = "phylo4d")}: show last few nodes} }
#' @author Ben Bolker, Thibaut Jombart
#' @seealso \code{\link{phylo4d}} constructor and \code{\linkS4class{phylo4d}}
#' class.
#' @keywords methods
setMethod("summary", signature(object="phylo4d"),
 function(object, quiet=FALSE) {
    x <- object
    res <- summary(as(x, "phylo4"), quiet=quiet)
    res$name <- deparse(substitute(object, sys.frame(-1)))
    tips <- tdata(object, "tip")
    nodes <- tdata(object, "internal")

    if (!quiet)
        cat("\nComparative data:\n")

    if (nrow(tips) > 0) {
        if(!quiet) {
            cat("\nTips: data.frame with", nTips(object), "taxa and",
                ncol(tips), "variable(s) \n\n")
        }
        sumry.tips <- summary(tips)
        res$sumry.tips <- sumry.tips
        if (!quiet)
            print(sumry.tips)
    }
    else {
        if (!quiet)
            cat("\nObject contains no tip data.")
    }
    if (nrow(nodes) > 0) {
        if (!quiet) {
            cat("\nNodes: data.frame with", nNodes(object), "internal nodes and",
                ncol(nodes), "variables \n\n")
        }
        sumry.nodes <- summary(nodes)
        res$sumry.nodes <- sumry.nodes
        if (!quiet)
            print(sumry.nodes)
    }
    else {
        if(!quiet)
            cat("\nObject contains no node data.\n")
    }
    invisible(res)
})


#' Tests for presence of data associated with trees stored as phylo4d objects
#' 
#' Methods that test for the presence of data associated with trees stored as
#' phylo4d objects.
#' 
#' The outcome of the test is based on row names of the data frame stored in
#' \code{data}. If there are no rows having row names from the set
#' \code{nodeId(x, "tip")}, then \code{hasTipData} returns FALSE.  Likewise, if
#' there are no rows having row names from the set \code{nodeId(x,
#' "internal")}, then \code{hasNodeData} returns FALSE.
#' 
#' @aliases hasNodeData hasNodeData-methods hasNodeData,phylo4d-method
#' hasTipData hasTipData-methods hasTipData,phylo4d-method
#' @param x a phylo4d object
#' @return \item{list("logical")}{return \code{TRUE} or \code{FALSE} depending
#' whether data are associated with the tree (i.e., the slots \code{tip.data}
#' or \code{node.data} are not empty)}
#' @section Methods: \describe{ \item{hasNodeData}{\code{signature(object =
#' "phylo4d")}: whether tree has internal node data}
#' \item{hasTipData}{\code{signature(object = "phylo4d")}: whether tree has
#' data associated with its tips} }
#' @author Ben Bolker, Thibault Jombart, Francois Michonneau
#' @seealso \code{\link{phylo4d}} constructor and \code{\linkS4class{phylo4d}}
#' class.
#' @keywords methods
#' @examples
#' 
#'   data(geospiza)
#'   hasTipData(geospiza)  ## TRUE
#'   hasNodeData(geospiza) ## FALSE
#' 
setMethod("hasTipData", signature(x="phylo4d"),
 function(x) {
    ncol(tdata(x, type="tip", empty.columns=FALSE)) > 0
})

setMethod("hasNodeData", signature(x="phylo4d"),
 function(x) {
    ncol(tdata(x, type="internal", empty.columns=FALSE)) > 0
})


## FIXME: doesn't deal with missing node data
##   (don't even know how that should be done in this case)
## setMethod("na.omit", "phylo4d", function(object, ...) {
##    tipdata <- tdata(object, "tip")
##    na.index <- which(!complete.cases(tipdata))
##    prune(object, tip = na.index)
##})
