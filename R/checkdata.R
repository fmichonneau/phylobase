
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
       length(object@label) == 0 && length(object@edge.label) == 0)
        return(TRUE)

    ## FIXME: check for cyclicity?
    nedges <- nrow(object@edge)

    if (hasEdgeLength(object)) {
      if (length(object@edge.length) != nedges)
        return("edge lengths do not match number of edges")
      if(!is.numeric(object@edge.length))
          return("edge lengths are not numeric")
      ## presumably we shouldn't allow NAs mixed
      ## with numeric branch lengths except at the root
      if (sum(is.na(object@edge.length)) > 1)
        return("NAs in edge lengths")
      ## Strip root edge branch length (if set to NA)
      if (any(object@edge.length[!is.na(object@edge.length)] < 0))
        return("edge lengths must be non-negative")
    }
    ##TODO fix this up somehow, or remove? (Nnode slot no longer exists)
    ## if (length(object@tip.label)+object@Nnode-1 != N) # does not work with multifurcations
    ##  return("number of tip labels not consistent with number of edges and nodes")
    ## check: tip numbers = (m+1):(m+n)
    ntips <- nTips(object)
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

    if (nRoots > 0) {
      if (sum(E[, 1] == 0) != 1) {
        return("for a rooted tree, edge matrix must contain (exactly one) explicit root edge with ancestor==0")
      }
      root.node <- unname(E[which(E[,1] == 0), 2])
      if (!root.node == nTips + 1)
        ## TODO this isn't actually a requirement
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

    ## make sure tip/node labels have internal names that match node IDs
    lab.msg <- "Use tipLabels<- (and nodeLabels<- if needed) to update them."
    if (is.null(names(object@label))) {
        return(c("Tip and node labels must have names matching node IDs. ",
            lab.msg))
             
    } else {
        if (!all(tips %in% names(na.omit(object@label)))) {
            return(c("All tips must have associated tip labels. ",
                lab.msg))
        }
        if (!all(names(object@label) %in% nodeId(object, "all"))) {
            return(c("One or more tip/node label has an unmatched ID name ",
                lab.msg))
        }
    }

    ## make sure edge lengths have internal names that match the edges
    elen.msg <- "Use edgeLength<- to update them."
    if(hasEdgeLength(object)) {
        if (is.null(names(object@edge.length))) {
            return(c("Edge lengths must have names matching edge IDs. ",
                elen.msg))
        }
        if (!all(names(object@edge.length) %in% edgeId(object, "all"))) {
            return(c("One or more edge lengths has an unmatched ID name. ",
                elen.msg))
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

# JR: I don't think this part is necessary. All that matters is that all
# rows in the data have names corresponding to (valid) node numbers
#    ntips <- nTips(object)
#    nnodes <- nNodes(object)
#
#    ## Check dimensions
#    if (nrow(object@tip.data) > 0 && nrow(object@tip.data) != ntips)
#        stop("The number of tip data does not match the number ",
#             "of tips in the tree")
#    if (nrow(object@node.data) > 0 && nrow(object@node.data) != nnodes)
#        stop("The number of node data does not match the number ",
#             "of internal nodes in the tree")

    ## Check rownames
    if (nrow(object@data) > 0 &&
        !all(row.names(object@data) %in% nodeId(object, "all")))
        stop("The row names of tree data do not match the node numbers")

    return(TRUE)
}
