## REQUIRED for all trees
checkPhylo4 <- function(object) {
    ct <- checkTree(object)

    if (class(object) == "phylo4d")
        ## checkPhyo4Data returns TRUE or fail
        cd <- checkPhylo4Data(object)

    return(ct)
}

checkTree <- function(object) {

    ## case of empty phylo4 object
    if(nrow(object@edge) == 0 && length(object@edge.length) == 0 &&
       length(object@label) == 0 && length(object@edge.label) == 0)
        return(TRUE)

    ## get options
    opt <- phylobase.options()

    ## Storage of error/warning messages
    err <- wrn <- character(0)

    ## Define variables
    nedges <- nEdges(object)
    ntips <- nTips(object)
    E <- edges(object)
    tips <- unique(sort(E[,2][!E[,2] %in% E[,1]]))
    nodes <- unique(sort(c(E)))
    intnodes <- nodes[!nodes %in% tips]
    roots <- E[which(is.na(E[,1])),2]
    nRoots <- length(roots)

    ## Check edge lengths
    if (hasEdgeLength(object)) {
        if (length(object@edge.length) != nedges)
            err <- c(err, "edge lengths do not match number of edges")
        if(!is.numeric(object@edge.length))
            err <- c(err, "edge lengths are not numeric")
        ## presumably we shouldn't allow NAs mixed
        ## with numeric branch lengths except at the root
        if (sum(is.na(object@edge.length)) > 1)
            err <- c(err, "NAs in edge lengths")
        ## Strip root edge branch length (if set to NA)
        if (any(object@edge.length[!is.na(object@edge.length)] < 0))
            err <- c(err, "edge lengths must be non-negative")
        ## Check edge length labels
        elen.msg <- "Use edgeLength<- to update them."
        if (is.null(names(object@edge.length))) {
            err <- c(err, paste("Edge lengths must have names matching edge IDs.",
                                elen.msg))
        }
        if (!all(names(object@edge.length) %in% edgeId(object, "all"))) {
            err <- c(err, paste("One or more edge lengths has an unmatched ID name.",
                                elen.msg))
        }
    }

    ## Make sure tips and
    if (!(all(tips==1:ntips) && all(nodes=(ntips+1):(ntips+length(intnodes)))))
        err <- c(err, "tips and nodes incorrectly numbered")

    ##careful - nAncest does not work for counting nRoots in unrooted trees
    nAncest <- tabulate(na.omit(E)[, 2],nbins=max(nodes)) ## bug fix from Jim Regetz
    nDesc <- tabulate(na.omit(E[,1]))
    nTips <- sum(nDesc==0)
    if (!all(nDesc[1:nTips]==0))
        err <- c(err, "nodes 1 to nTips must all be tips")

    if (nRoots > 0) {
      if (sum(E[, 1] == 0) != 1) {
          err <- c(err, "for a rooted tree, edge matrix must contain (exactly one) explicit root edge with ancestor==0")
      }
      root.node <- unname(E[which(E[,1] == 0), 2])
    }

    ## Check that nodes are correctly numbered
    if (!all(nDesc[(nTips+1):(nTips+nNodes(object))]>0))
        err <- c(err, "nodes (nTips+1) to (nTips+nNodes) must all be internal nodes")

    ## how do we identify loops???
    ## EXPERIMENTAL: could be time-consuming for large trees?
    if (FALSE) {
      Emat <- matrix(0,nrow=max(E),ncol=max(E))
      Emat[E] <- 1
    }
    if (!object@order %in% phylo4_orderings) {
      err <- c(err, paste("unknown order: allowed values are",
               paste(phylo4_orderings,collapse=",")))
    }

    ## make sure tip/node labels have internal names that match node IDs
    lab.msg <- "Use tipLabels<- (and nodeLabels<- if needed) to update them."
    if (is.null(names(object@label))) {
        err <- c(err, paste("Tip and node labels must have names matching node IDs.",
                            lab.msg))

    } else {
        if (!all(tips %in% names(na.omit(object@label)))) {
            err <- c(err, paste("All tips must have associated tip labels.",
                                lab.msg))
        }
        if (!all(names(object@label) %in% nodeId(object, "all"))) {
            err <- c(err, paste("One or more tip/node label has an unmatched ID name",
                                lab.msg))
        }
    }

    ## make sure edge labels have internal names that match the edges
    elab.msg <- "Use edgeLabels<- to update them."
    if(hasEdgeLabels(object)) {
        if (is.null(names(object@edge.label))) {
            err <- c(err, paste("Edge labels must have names matching edge IDs.",
                                elab.msg))
        }
        if (!all(names(object@edge.label) %in% edgeId(object, "all"))) {
            err <- c(err, paste("One or more edge labels has an unmatched ID name.",
                                elab.msg))
        }
    }

    ## make sure that tip and node labels are unique
    if (hasDuplicatedLabels(object)) {
        currmsg <- "Labels are not unique"
        if (opt$allow.duplicated.labels == "fail")
            err <- c(err, currmsg)
        if (opt$allow.duplicated.labels == "warn")
            wrn <- c(wrn, currmsg)
    }

    if (any(nDesc>2)) {
        currmsg <- "tree includes polytomies"
        if (opt$poly == "fail")
            err <- c(err, currmsg)
        if (opt$poly == "warn")
            wrn <- c(wrn, currmsg)
      }

    if (nRoots>1) {
        currmsg <- "tree has more than one root"
        if (opt$multiroot == "fail")
            err <- c(err, currmsg)
        if (opt$multiroot == "warn")
            wrn <- c(wrn,currmsg)
    }
    if (any(nDesc==1)) {
        currmsg <- "tree contains singleton nodes"
        if (opt$singleton == "fail")
            err <- c(err, currmsg)
        if (opt$singleton == "warn")
            wrn <- c(wrn, currmsg)
    }
    if (any(nAncest>1)) {
      currmsg <- paste("tree is reticulated [most functions in phylobase haven't",
                       "been tested with reticulated trees]")
      if (opt$retic == "fail")
          err <- c(err, currmsg)
      if (opt$retic == "warn")
          wrn <- c(wrn, currmsg)
    }
    if (length(wrn) > 0) {
        wrn <- paste(wrn, collapse=", ")
        warning(wrn)
    }
    if (length(err) > 0) {
        err <- paste(err, collapse=", ")
        return(err) #failures are returned as text
    }
    else {
        return(TRUE)
    }
}

checkPhylo4Data <- function(object) {

    ## These are just some basic tests to make sure that the user does not
    ## alter the object in a significant way

    ## Check rownames
    if (nrow(object@data) > 0 &&
        !all(row.names(object@data) %in% nodeId(object, "all")))
        stop("The row names of tree data do not match the node numbers")

    return(TRUE)
}
