
## REQUIRED for all trees
checkPhylo4 <- function(object) {
    ct <- checkTree(object)

    if (class(object) == "phylo4d")
        ## checkPhyo4Data returns TRUE or fail
        cd <- checkPhylo4Data(object)

    return(ct)
}

checkTree <- function(object,
                      warn=c("retic","singleton","multiroot"),
                      err=NULL) {

    ## case of empty phylo4 object
    if(nrow(object@edge) == 0 && length(object@edge.length) == 0 &&
       object@Nnode == 0 && length(object@node.label) == 0 &&
       length(object@tip.label) == 0 && length(object@edge.label) == 0)
        return(TRUE)

    ## FIXME: check for cyclicity?
    nedges <- nrow(object@edge)

    if (hasEdgeLength(object)) {
      if (length(object@edge.length) != nedges)
        return("edge lengths do not match number of edges")
      if(!is.numeric(object@edge.length))
          stop("Edge lengths are not numeric.")
      ## presumably we shouldn't allow NAs mixed
      ## with numeric branch lengths except at the root
      if (sum(is.na(object@edge.length)) > 1)
        return("NAs in edge lenghts")
      ## Strip root edge branch lenght (if set to NA)
      if (any(object@edge.length[!is.na(object@edge.length)] < 0))
        return("edge lengths must be non-negative")
    }
    ## if (length(object@tip.label)+object@Nnode-1 != N) # does not work with multifurcations
    ##  return("number of tip labels not consistent with number of edges and nodes")
    ## check: tip numbers = (m+1):(m+n)
    ntips <- nTips(object)
    if(length(object@tip.label) != ntips)
      return("number of tip labels not consistent with number of tips")
    E <- edges(object)
    tips <- unique(sort(E[,2][!E[,2] %in% E[,1]]))
    nodes <- unique(sort(c(E)))
    intnodes <- nodes[!nodes %in% tips]
    roots <- E[which(is.na(E[,1])),2]
    nRoots <- length(roots)

    if (!(all(tips==1:ntips) && all(nodes=(ntips+1):(ntips+length(intnodes)))))
      return("tips and nodes incorrectly numbered")

    ##careful - nAncest does not work for counting nRoots in unrooted trees
    nAncest <- tabulate(na.omit(E)[, 2],nbins=max(nodes)) ## bug fix from Jim Regetz
    nDesc <- tabulate(na.omit(E[,1]))
    nTips <- sum(nDesc==0)
    if (!all(nDesc[1:nTips]==0))
      return("nodes 1 to nTips must all be tips")

    if (nRoots>0) {
      if (sum(is.na(E[,1]))!=1) {
        return("for a rooted tree, edge matrix must contain (exactly one) explicit root edge with ancestor==NA")
      }
      root.node <- unname(E[which(is.na(E[,1])),2])
      if (!root.node==nTips+1)
        return("root node must be first row of edge matrix")
    }

    ##
    ## how do we identify loops???
    ## EXPERIMENTAL: could be time-consuming for large trees?
    if (FALSE) {
      Emat <- matrix(0,nrow=max(E),ncol=max(E))
      Emat[E] <- 1
    }
    if (!object@order %in% phylo4_orderings) {
      stop("unknown order: allowed values are ",
           paste(phylo4_orderings,collapse=","))
    }

    ## make sure that nodes and edges have internal names
    ## and that they match the nodes
    if (is.null(names(object@tip.label))) {
        if(length(object@tip.label) == nTips(object)) {
            stop("There is no internal name associated with your tips. Use the ",
                 "function tipLabels <- to change your tip labels.")
        }
        else
            stop("Your object doesn't have internal node names and the number of ",
                 "tip labels doesn't match the number tips.")
    }
    else {
        if(!all(names(object@tip.label) %in%  nodeId(object, "tip")))
            stop("Internal names for tips don't match tip ID numbers")
    }

    if (is.null(names(object@node.label))) {
        if(length(object@node.label) == nNodes(object)) {
            stop("There is no internal names associated with internal ",
                 "nodes. Use the function nodeLabels <- to create or ",
                 "change your internal node labels.")
        }
        else
            stop("Your object doesn't have internal node names and the number of ",
                 "node labels doesn't match the number nodes.")
    }
    else {
        if(!all(names(object@node.label) %in%  nodeId(object, "internal")))
            stop("Internal names for tips don't match tip ID numbers")
    }

    if(hasEdgeLength(object)) {
        if(is.null(names(object@edge.length))) {
            warning("Your edges don't have internal names. Use the function ",
                    "edgeLength <- to update the the branch lengths of your ",
                    "tree.")
        }
        else {
            tEdgLbl <- paste(object@edge[,1], object@edge[,2], sep="-")
            if(!all(names(object@edge.length) %in% tEdgLbl))
                stop("There is something wrong with your internal edge length ",
                     "labels. Use the function edgeLength <- to update the the ",
                     "branch lengths of your tree.")
        }
    }

    ## make sure that tip and node labels are unique
    lb <- labels(object, "allnode")
    lb <- lb[nchar(lb) > 0]
    lb <- na.omit(lb)
    if(any(table(lb) > 1))
        stop("All labels must be unique")

    ## all done with fatal errors.  Now construct a list
    ##  of warnings and paste them together
    msg <- character(0)

    ##fixme following check fails for unrooted trees
    ##if (!all(nDesc[(nTips+1):(nTips+nNodes(object))]>0))
    ##  return("nodes (nTips+1) to (nTips+nNodes) must all be internal nodes")
    if (any(nDesc>2)) {
        currmsg <- "tree includes polytomies"
        if ("poly" %in% err)
          return(currmsg)
        if ("poly" %in% warn)
          msg <- c(msg,currmsg)
      }

    if (nRoots>1) {
        currmsg <- "tree has more than one root"
        if ("multiroot" %in% err)
          return(currmsg)
        if ("multiroot" %in% warn)
          msg <- c(msg,currmsg)
      }
    if (any(nDesc==1)) {
        currmsg <- "tree contains singleton nodes"
          if ("singleton" %in% err)
            return(currmsg)
          if ("singleton" %in% warn)
            msg <- c(msg,currmsg)
      }
    if (any(nAncest>1)) {
      currmsg <- paste("tree is reticulated [most functions in phylobase haven't",
                       "been tested with reticulated trees]")
      if ("retic" %in% err)
        return(currmsg)
      if ("retic" %in% warn)
        msg <- c(msg,currmsg)
    }
    if (length(msg)>0) {
      msg <- paste(msg,collapse=", ")
      warning(msg)
    }
    return(TRUE)
  }

checkPhylo4Data <- function(phy) {

    ## These are just some basic tests to make sure that the user does not
    ## alter the object in a significant way

    ntips <- nTips(phy)
    nnodes <- nNodes(phy)

    ## Check dimensions
    if (nrow(phy@tip.data) > 0 && nrow(phy@tip.data) != ntips)
        stop("The number of tip data does not match the number ",
             "of tips in the tree")
    if (nrow(phy@node.data) > 0 && nrow(phy@node.data) != nnodes)
        stop("The number of node data does not match the number ",
             "of internal nodes in the tree")

    ## Check rownames
    if (nrow(phy@tip.data) > 0 &&
       !all(rownames(phy@tip.data) %in% nodeId(phy, "tip")))
        stop("The row names of tip data do not match the tip numbers")
    if (nrow(phy@node.data) > 0 &&
        !all(rownames(phy@node.data) %in% nodeId(phy, "internal")))
        stop("The row names of node data do not match the node numbers")

    return(TRUE)
}



formatData <- function(phy, dt, type=c("tip", "internal", "all"),
                       match.data=TRUE, rownamesAsLabels=FALSE,
                       label.type=c("rownames", "column"),
                       label.column=1, missing.data=c("fail", "warn", "OK"),
                       extra.data=c("warn", "OK", "fail")
                       ) {

    type <- match.arg(type)
    label.type <- match.arg(label.type)
    stopifnot(label.column %in% 1:ncol(dt))
    missing.data <- match.arg(missing.data)
    extra.data <- match.arg(extra.data)

    nr <- switch(type,
                 tip = nTips(phy),
                 internal = nNodes(phy),
                 all = nTips(phy)+nNodes(phy))

    tmpDt <- array(, dim=c(nr, ncol(dt)),
                   dimnames=list(nodeId(phy, type), colnames(dt)))
    tmpDt <- data.frame(tmpDt)

    if(match.data) {
        ## Replace node labels by node numbers
        ndNames <- switch(label.type,
                          rownames = rownames(dt),
                          column = dt[,label.column])
        ndDt <- lapply(ndNames, function(nd) {
            if(nchar(gsub("[0-9]", "", nd)) == 0 && !rownamesAsLabels)
                getNode(phy, as.integer(nd), missing="OK")
            else getNode(phy, nd, missing="OK")
        })
        ndDt <- unlist(ndDt)

        ## Make sure that data are matched to appropriate nodes
        if(type != "all") {
            switch(type,
                   tip = {
                     ## BMB: don't bother trying to match NAs
                       if(any(na.omit(names(ndDt)) %in% labels(phy, "internal")))
                           stop("You are trying to match tip data to internal ",
                                "nodes. Make sure that your data identifiers ",
                                "are correct.")
                   },
                   internal = {
                       if(any(na.omit(names(ndDt)) %in% labels(phy, "tip")))
                           stop("You are trying to match node data to tip ",
                                "nodes. Make sure that your data identifiers ",
                                "are correct.")
                   })
        }

        ## Check differences
        extra <- names(ndDt[is.na(ndDt)])
        mssng <- nodeId(phy, type)[! nodeId(phy, type) %in% ndDt]

        if(length(mssng) > 0 && missing.data != "OK") {
            msg <- "The following nodes are not found in the dataset: "

            ## provides label if it exists and node number otherwise
            mssng <- sapply(mssng, function(m) {
                m <- getNode(phy, m)
                if (is.na(names(m)) || is.null(names(m)))
                    m
                else
                    names(m)
            })

            msg <- paste(msg, paste(mssng, collapse=", "))
            switch(missing.data,
                   warn = warning(msg),
                   fail = stop(msg))
        }

        if(length(extra) > 0 && extra.data != "OK") {
            msg <- "The following names are not found in the tree: "
            msg <- paste(msg, paste(extra, collapse=", "))
            switch(extra.data,
                   warn = warning(msg),
                   fail = stop(msg))

        }
        ## Format data to have correct dimensions
        dt <- dt[!is.na(ndDt) ,, drop=FALSE]
        rownames(dt) <- ndDt[!is.na(ndDt)]
        if(label.type == "column") dt <- dt[, -label.column]
        tmpDt <- dt[match(rownames(tmpDt), rownames(dt)) ,, drop=FALSE]
    }
    else {
        ## Remove rownames in data provided
        rownames(dt) <- NULL

        ## Check differences between dataset and tree
        diffNr <- nrow(dt) - nr
        if(diffNr > 0 && extra.data != "OK") {
            msg <- paste("There are", diffNr, "extra rows.")
            switch(extra.data,
                   warn = warning(msg),
                   fail = stop(msg))
        }
        if(diffNr < 0 && missing.data != "OK") {
            msg <- paste("There are", abs(diffNr), "missing rows.")
            switch(missing.data,
                   warn = warning(msg),
                   fail = stop(msg))
        }
        tmpDt <- dt[1:min(nrow(dt), nr) ,, drop = FALSE]
    }

    tmpDt
}
