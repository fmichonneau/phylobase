setMethod("print", "phylo4d", printphylo4)

setMethod("show", "phylo4d", function(object) printphylo4(object))

setMethod("tdata", "phylo4d", function(x, which = c("tip",
    "node", "allnode"), label.type=c("row.names","column"), ...) {
    which <- match.arg(which)
    ## FIXME: should have no labels in this case?
    if (!hasNodeLabels(x) && which=="node" && missing(label.type)) { }
    label.type <- match.arg(label.type)
    if (which == "tip") {
        if (all(dim(x@tip.data)==0)) {
            return(x@tip.data)
        }
        tdata <- x@tip.data
        data.names <- x@tip.label
        if ( identical(label.type,"row.names") ) {
            if ( identical(data.names,unique(data.names)) ||
                !(any(is.na(data.names))) ) {
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
            tdata <- data.frame(label=data.names,tdata)
        }
        return(tdata)
    }

    if (which == "node") {
        if (all(dim(x@node.data)==0)) {
            return(x@node.data)
        }
        tdata <- x@node.data
        data.names <- x@node.label
        if ( identical(label.type,"row.names") ) {
          if ( length(data.names)>0 &&
              !any(duplicated(data.names)) &&
              !(any(is.na(data.names)))) {
            row.names(tdata) <- data.names
          } else {
            warning("Non-unique or missing labels found, ",
                    "labels cannot be coerced to tdata row.names. ",
                    "Use the label.type argument to include labels ",
                    "as first column of data.")
          }
        }
        if (identical(label.type,"column")) {
          if (!hasNodeLabels(x)) data.names <- rep("",nNodes(x))
          tdata <- data.frame(label=data.names,tdata)
        }
        return(tdata)
    }

    if (which == "allnode") {
        if (all(dim(x@node.data)==0)) { ## empty data
          if (!hasNodeLabels(x)) {
            nodedata <- data.frame(label=rep("",nNodes(x)))
          } else
          nodedata <- data.frame(label=x@node.label)
        } 
        else {
          nodedata <- tdata(x, "node", label.type="column")
        }
        if (all(dim(x@tip.data)==0)) {
          tipdata <- data.frame(label=x@tip.label)
        }
        else {
            tipdata <- tdata(x, "tip", label.type="column")
        }

        data.names <- c(as.character(nodedata$label),as.character(tipdata$label))
        tipdata$label <- (x@Nnode+1):(x@Nnode+length(x@tip.label))
        nodedata$label <- 1:x@Nnode
        ## FIXME - kludgy merge and subsequent cleanup - make robust
        tdata <- merge(nodedata,tipdata, all=TRUE,sort=FALSE)[,-1,drop=FALSE]
        tdata <- data.frame(label=data.names,tdata)

        if ( identical(label.type,"row.names") ) {
            if ( identical(data.names,unique(data.names)) || !(any(is.na(data.names))) ) {
                tdata <- data.frame(tdata[,-1,drop=FALSE])
                row.names(tdata) <- data.names
            }
            else {
                stop("Non-unique or missing labels found, labels cannot be coerced to tdata row.names. Use the label.type argument to include labels as first column of data.")
            }
        }
        return(tdata)
    }
})

setReplaceMethod("tdata", "phylo4d", function(object, which = c("tip",
    "node", "allnode"), ..., value) {
    which <- match.arg(which)
    if (which == "allnode") {
        namesmatch <- all(colnames(object@tip.data) == colnames(object@node.data))
        classmatch <- all(sapply(object@tip.data, class) == sapply(object@node.data,
            class))
        if (!(classmatch && namesmatch))
            stop("Node and tip columns do not match;",
                 "you should access tip and node data separately")
    }
    if(is.matrix(value)) value <- as.data.frame(value)
    if(!is.data.frame(value))
        stop("For now, only data.frame or matrix can be provided")
    switch(which,
           tip = object@tip.data <- value,
           node = object@node.data <- value,
           allnode = stop("for now, must set tip and node data separately"))
    if(check_data(object, ...)) object <- attach_data(object, ...)
    object
})


## Alternative phylo4d summary method, using phylo4 summary
## Marguerite Butler & Peter Cowan
setMethod("summary", "phylo4d", function(object) {
    x <- object
    summary(as(object, "phylo4"))
    tips <- tdata(object, "tip")
    nodes <- tdata(object, "node")
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

setReplaceMethod("nodeLabels", "phylo4d", function(object, ...,
    value) {
    object@node.label <- value
    #rownames(object@node.data) <- value
    object
})

setReplaceMethod("labels", "phylo4d", function(object, ..., value) {
    object@tip.label <- value
    #rownames(object@tip.data) <- value
    object
})

## FIXME: doesn't deal with missing node data
##   (don't even know how that should be done in this case)
setMethod("na.omit", "phylo4d", function(object, ...) {
    tipdata <- tdata(object, "tip")
    na.index <- which(!complete.cases(tipdata))
    prune(object, tip = na.index)
})

setMethod("names", signature(x = "phylo4d"), function(x) {
    temp <- rev(names(attributes(x)))[-1]
    return(rev(temp))
})

setMethod("reorder", signature(x = 'phylo4d'), function(x, order = 'cladewise') {
        index <- orderIndex(x, order)
        test <<- index
        x@edge      <- x@edge[index, ]
        x@tip.data  <- x@tip.data[index[index <= nTips(x)], , drop = FALSE]
        x@node.data <- x@node.data[index[index > nTips(x)], , drop = FALSE]
        if(hasEdgeLabels(x)) { x@edge.label  <- x@edge.label[index] }
        if(hasEdgeLength(x)) { x@edge.length <- x@edge.length[index] }
        x
    })
