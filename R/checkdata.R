
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
        return("NAs in edge lengths")
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
    lb <- labels(object, "all")
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

checkPhylo4Data <- function(object) {

    ## These are just some basic tests to make sure that the user does not
    ## alter the object in a significant way

    ntips <- nTips(object)
    nnodes <- nNodes(object)

    ## Check dimensions
    if (nrow(object@tip.data) > 0 && nrow(object@tip.data) != ntips)
        stop("The number of tip data does not match the number ",
             "of tips in the tree")
    if (nrow(object@node.data) > 0 && nrow(object@node.data) != nnodes)
        stop("The number of node data does not match the number ",
             "of internal nodes in the tree")

    ## Check rownames
    if (nrow(object@tip.data) > 0 &&
       !all(rownames(object@tip.data) %in% nodeId(object, "tip")))
        stop("The row names of tip data do not match the tip numbers")
    if (nrow(object@node.data) > 0 &&
        !all(rownames(object@node.data) %in% nodeId(object, "internal")))
        stop("The row names of node data do not match the node numbers")

    return(TRUE)
}
