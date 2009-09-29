setMethod("tdata", signature(x="phylo4d"),
  function(x, type=c("tip", "internal", "all"),
           label.type=c("row.names","column"),
           empty.columns=TRUE, ...) {

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
 function(x, type = c("tip", "internal", "all"), ..., value) {
    type <- match.arg(type)
    object <- x

    ## Removes existing data, just keeps the tree (as a phylo4d)
    object <- extractTree(object)
    object <- as(object, "phylo4d")

    object@data <- switch(type,
                      tip = .phylo4Data(object, tip.data=value, ...),
                      internal = .phylo4Data(object, node.data=value, ...),
                      all = .phylo4Data(object, all.data=value, ...))

    object
})

setMethod("addData", signature(x="phylo4d"),
  function(x, tip.data=NULL, node.data=NULL,
           all.data=NULL, pos=c("after", "before"),
           merge.data=TRUE, match.data=TRUE,  ...) {

    pos <- match.arg(pos)

    new.data <- .phylo4Data(x=x, tip.data=tip.data, node.data=node.data,
        all.data=all.data, merge.data=merge.data, match.data=match.data, ...)

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
  function(x, tip.data=NULL, node.data=NULL,
           all.data=NULL, pos=c("after", "before"),
           merge.data=TRUE, match.data=TRUE, ...) {
    phylo4d(x, tip.data=tip.data, node.data=node.data, all.data=all.data,
            merge.data=merge.data, match.data=match.data, ...)
})

## Alternative phylo4d summary method, using phylo4 summary
## Marguerite Butler & Peter Cowan
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

#setMethod("tipData", signature(x="phylo4d"),
# function(x) {
#    nrow(x@tip.data) > 0
#})

setMethod("hasTipData", signature(x="phylo4d"),
 function(x) {
    ncol(tdata(x, type="tip", empty.columns=FALSE)) > 0
})

#setMethod("nodeData", signature(x="phylo4d"),
# function(x) {
#    nrow(x@tip.data) > 0
#})

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
