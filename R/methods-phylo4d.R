setMethod("tdata", "phylo4d",
  function(x, type=c("tip", "internal", "allnode"),
           label.type=c("row.names","column"),
           empty.columns=TRUE, ...) {

      ## Returns data associated with the tree
      ## Note: the function checks for unique labels. It's currently unecessary
      ## but could be useful in the future if non-unique labels are allowed.

      type <- match.arg(type)
      label.type <- match.arg(label.type)

      if (type == "tip") {
          if (all(dim(x@tip.data) == 0)) {
              return(x@tip.data)
          }
          tdata <- x@tip.data
          data.names <- tipLabels(x)
          if ( label.type ==  "row.names" ) {
              if (!any(duplicated(data.names)) &&
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
      }

      if (type == "internal") {
          if (all(dim(x@node.data)==0)) {
              return(x@node.data)
          }
          tdata <- x@node.data
          if(hasNodeLabels(x))
              data.names <- nodeLabels(x)
          else
              data.names <- nodeId(x, "internal")

          if ( identical(label.type, "row.names") ) {
              if ( length(data.names) > 0 &&
                  !any(duplicated(data.names)) &&
                  !(any(is.na(data.names)))) {
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
      }

      if (type == "allnode") {
          ## node data
          if (all(dim(x@node.data) == 0)) { # empty data
              if (!hasNodeLabels(x)) {
                  nodedata <- data.frame(label=x@node.label)
              }
              else
                  nodedata <- data.frame(label=nodeLabels(x))
          }
          else {
              nodedata <- tdata(x, "internal", label.type="column")
          }

          ## tip data
          if (all(dim(x@tip.data) == 0)) {
              tipdata <- data.frame(label=tipLabels(x))
          }
          else {
              tipdata <- tdata(x, "tip", label.type="column")
          }

          ## following lines necessary to be able to use merge on data
          ## belonging to different classes (e.g. nodeId as numeric and
          ## labels as character)
          tipdata$label <- as.character(tipdata$label)
          nodedata$label <- as.character(nodedata$label)

          tdata <- merge(tipdata, nodedata, all=TRUE, sort=FALSE)[,, drop=FALSE]

          if (identical(label.type, "row.names")) {
              if (identical(tdata$label, unique(tdata$label)) ||
                  !(any(is.na(tdata$label))) ) {
                  row.names(tdata) <- tdata[,1]
                  tdata <- data.frame(tdata[, -1, drop=FALSE])
              }
              else {
                  stop("Non-unique or missing labels found, labels cannot be ",
                       "coerced to tdata row.names. Use the label.type argument ",
                       "to include labels as first column of data.")
              }
          }

      }

      ## remove empty columns (filled with NAs)
      if(!empty.columns) {
          emptyCol <- apply(tdata, 2, function(x) all(is.na(x)))
          tdata <- tdata[, !emptyCol]
      }

      tdata
  })

setReplaceMethod("tdata", "phylo4d",
                 function(object, type = c("tip", "internal", "allnode"), ...,
                          value) {
    type <- match.arg(type)

    tmpData <- switch(type,
                      tip = .phylo4Data(object, tip.data=value, ...),
                      internal = .phylo4Data(object, node.data=value, ...),
                      allnode = .phylo4Data(object, all.data=value, ...))

    if(all(dim(tmpData$tip.data)))
        object@tip.data <- tmpData$tip.data
    if(all(dim(tmpData$node.data)))
        object@node.data <- tmpData$node.data

    object
})

setMethod("addData", "phylo4d", function(x, tip.data=NULL, node.data=NULL,
                                         all.data=NULL, pos=c("after", "before"),
                                         merge.data=TRUE, match.data=TRUE,
                                         ...) {
    pos <- match.arg(pos)

    tmpData <- .phylo4Data(x=x, tip.data=tip.data, node.data=node.data,
                           all.data=all.data, merge.data=merge.data,
                           match.data=match.data, ...)

    if(identical(pos, "before")) {
        if(!all(dim(tmpData$tip.data) == 0)) {
            if(all(dim(x@tip.data) > 0))
                x@tip.data <- cbind(tmpData$tip.data, x@tip.data)
            else
                x@tip.data <- tmpData$tip.data
        }
        if(!all(dim(tmpData$node.data) == 0)) {
            if(all(dim(x@tip.data) > 0))
                x@node.data <- cbind(tmpData$node.data, x@node.data)
            else
                x@node.data <- tmpData$node.data
        }
    }
    else {
        if(!all(dim(tmpData$tip.data) == 0)) {
            if(all(dim(x@tip.data) > 0))
                x@tip.data <- cbind(x@tip.data, tmpData$tip.data)
            else
                x@tip.data <- tmpData$tip.data
        }

        if(!all(dim(tmpData$node.data) == 0)) {
            if(all(dim(x@node.data) > 0))
                x@node.data <- cbind(x@node.data, tmpData$node.data)
            else
                x@node.data <- tmpData$node.data
        }
    }

    x
})

setMethod("addData", "phylo4", function(x, tip.data=NULL, node.data=NULL,
                                        all.data=NULL, pos=c("after", "before"),
                                        merge.data=TRUE, match.data=TRUE,
                                        ...) {
    phylo4d(x, tip.data=tip.data, node.data=node.data, all.data=all.data,
            merge.data=merge.data, match.data=match.data, ...)
})

## Alternative phylo4d summary method, using phylo4 summary
## Marguerite Butler & Peter Cowan
setMethod("summary", "phylo4d", function(object) {
    x <- object
    summary(as(x, "phylo4"))
    tips <- tdata(object, "tip")
    nodes <- tdata(object, "internal")
    cat("\nComparative data:\n")
    if (nrow(tips) > 0) {
        cat("\nTips: data.frame with", nTips(object), "taxa and",
            ncol(tips), "variable(s) \n\n")
        print(summary(tips))
    }
    else {
        cat("\nObject contains no tip data.")
    }
    if (nrow(nodes) > 0) {
        cat("\nNodes: data.frame with", nNodes(object), "internal nodes and",
            ncol(nodes), "variables \n\n")
        print(summary(nodes))
    }
    else {
        cat("\nObject contains no node data.\n")
    }
})

setMethod("hasNodeData", "phylo4d", function(x) {
    nrow(x@node.data) > 0
})


## FIXME: doesn't deal with missing node data
##   (don't even know how that should be done in this case)
setMethod("na.omit", "phylo4d", function(object, ...) {
    tipdata <- tdata(object, "tip")
    na.index <- which(!complete.cases(tipdata))
    prune(object, tip = na.index)
})
