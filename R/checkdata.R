
## REQUIRED for all trees
check_phylo4 <- function(object) {
    check_tree(object)
}

check_tree <- function(object,warn="retic",err=NULL) {
    ## FIXME: check for cyclicity?
    N <- nrow(object@edge)  
    if (hasEdgeLength(object) && length(object@edge.length) != N)
      return("edge lengths do not match number of edges")
    ## if (length(object@tip.label)+object@Nnode-1 != N) # does not work with multifurcations
    ##  return("number of tip labels not consistent with number of edges and nodes")
    ## check: internal node numbers = 1:m
    
    ## check: tip numbers = (m+1):(m+n)
    ntips <- nTips(object)
    if(length(object@tip.label) != ntips)
      return("number of tip labels not consistent with number of tips")
    E <- edges(object)
    tips <- sort(E[,2][!E[,2] %in% E[,1]])
    nodes <- unique(sort(c(E)))
    intnodes <- nodes[!nodes %in% tips]
    if (!(all(tips==1:ntips) && all(nodes=(ntips+1):(ntips+length(intnodes)))))
      return("tips and nodes incorrectly numbered")
    nAncest <- tabulate(E[, 2])
    nDesc <- tabulate(E[,1])
    nTips <- sum(nDesc==0)
    if (!all(nDesc[1:nTips]==0))
      return("nodes 1 to nTips must all be tips")
    if (!all(nDesc[(nTips+1):(nTips+nNodes(object))]>0))
      return("nodes (nTips+1) to (nTips+nNodes) must all be internal nodes")
    if (any(nDesc>2)) {
        if ("poly" %in% err)
          return("tree includes polytomies")
        if ("poly" %in% warn)
          warning("tree includes polytomies")
    }
    nRoots <- sum(nAncest==0)
    if (which(nAncest==0)!=nTips+1) {
      return("root node is not at position (nTips+1)")
    }
    if (any(nAncest==0) && E[1,1]!=nTips+1) {
      return("root node must be first row of edge matrix")
    }
    ##
    ## how do we identify loops???
    ## EXPERIMENTAL: could be time-consuming for large trees?
    if (FALSE) {
      Emat <- matrix(0,nrow=max(E),ncol=max(E))
      Emat[E] <- 1
    }
    ## all done with fatal errors.  Now construct a list
    ##  of warnings and paste them together
    msg <- character(0)
    if (nRoots>1)
      msg <- "tree has more than one root"
    ## BMB: should this be an error????
    if (any(nAncest>1))
      msg <- c(msg,"some nodes have multiple ancestors")
    if (any(nDesc==1))
      msg <- c("tree contains singleton nodes")
    msg <- paste(msg,collapse=", ")
    if (nzchar(msg)) {
        if ("retic" %in% err)
          return(paste("tree is reticulated:",msg))
        if ("retic" %in% warn)
          warning("tree is reticulated:",msg)
    }
    return(TRUE)
}


check_data <- function(object,
                       label.type=c("row.names","column"),
                       label.column=1,
                       use.tip.names=TRUE,
                       missing.tip.data=c("fail","OK","warn"),
                       extra.tip.data=c("fail","OK","warn"),
                       default.tip.names=c("warn","OK","fail"),
                       use.node.names=FALSE,
                       missing.node.data=c("OK","warn","fail"),		
                       extra.node.data=c("OK","warn","fail"),												
                       default.node.names=c("warn","OK","fail"),...)							 

{

    ## name matching default: use row.names of data frame
    label.type = match.arg(label.type)
    if (identical(label.type, "row.names")) {
        tip.names <- row.names(object@tip.data)
        node.names <- row.names(object@node.data)
    }
    else {
        tip.names <- object@tip.data[,label.column]
        node.names <- object@node.data[,label.column]        
    }
    
    ## tip default: use names, require names, must match exactly
    missing.tip.data <- match.arg(missing.tip.data)
    extra.tip.data <- match.arg(extra.tip.data)
    default.tip.names <- match.arg(default.tip.names)
    
    ## node default: don't use node names, don't require names, do not need to match exactly
    missing.node.data <- match.arg(missing.node.data)
    extra.node.data <- match.arg(extra.node.data)
    default.node.names <- match.arg(default.node.names)
    
    ## for each set of data, check for names, missing and extra data and take appropriate actions
    
    ## tip data checks
    ## if tip.data exist
    if (!all(dim(object@tip.data)==0)) {
        ## if we want to use tip.names
        if (use.tip.names) {
            
            ## check for default names
            if (all(tip.names == 1:length(tip.names))) {
                ## no tip.names
                if (default.tip.names == "fail") {
                    stop("Tip data have default names and may not match tree tip labels. ",
                         "Consider using the use.tip.names=FALSE option.")
                }
                else if (default.tip.names == "warn") {
                    warning("Tip data have default names and may not match tree tip labels. ",
                            "Consider using the use.tip.names=FALSE option.")
                }
            }
            
            ## check tip names
            ## check for missing or extra tip data (relative to tree taxa)
            if (setequal(tip.names, object@tip.label)) {
                ## names are perfect match - ok
                return(TRUE)
            }
            else {
                ## we know the tree taxa and tip.data taxa are not a perfect match
                ## if tip.data taxa are subset of tree taxa, check missing.tip.data arg and act accordingly
                tips.in.rownames <- object@tip.label %in% tip.names
                rownames.in.tips <- tip.names %in% object@tip.label 
                missing.data.names <- object@tip.label[!tips.in.rownames]
                missing.data.name.msg <- if (length(missing.data.names)==0) "" else {
                    paste("\n(missing data names: ",
                          paste(missing.data.names,collapse=","),")",sep="")
                }
                extra.data.names <- tip.names[!rownames.in.tips]
                extra.data.name.msg <- if (length(extra.data.names)==0) "" else {
                    paste("\n(extra data names: ",
                          paste(extra.data.names,collapse=","),")",sep="")
                }
                if (!all(tips.in.rownames)) {
                    ## we know it's not an exact match - we have missing.tip.data - take action
                    if (!any(tips.in.rownames)) {
                        errmsg <- paste("No matches between tip data names and tree tip labels.",
                                        missing.data.name.msg,extra.data.name.msg)
                        if (missing.tip.data == "fail") {
                            stop(errmsg)
                        }
                        else if (missing.tip.data == "warn") {
                            warning(errmsg)
                        }
                    }
                    else
                      {
                          errmsg <- paste("Tip data names are a subset of tree tip labels",
                                          missing.data.name.msg,
                                          extra.data.name.msg)
                          if (missing.tip.data == "fail") {
                              stop(errmsg)
                          }

                          else if (missing.tip.data == "warn") {
                              warning(errmsg)
                          }
                      }
                    ##else ok
                }
                
                ##if tree taxa are subset of tip.data, check extra.tip arg and act accordingly
                if (!all(tip.names %in% object@tip.label)) {
                    ##we know it's not an exact match - we have extra.tip.data - take action
                    ##fail
                    errmsg <- paste("Tip data names are a superset of tree tip labels",
                                    missing.data.name.msg,
                                    extra.data.name.msg)
                    if (extra.tip.data == "fail") {
                        stop(errmsg)
                    }
                    ##warn
                    else if (extra.tip.data == "warn") {
                        warning(errmsg)
                    }
                    ##else ok
                }
                
                return(TRUE)
            } 
        }
        else
          {
              ##don't use tip names or attempt to sort - but check to make sure dimensions match
              if (!(nTips(object)==dim(object@tip.data)[1])) {
                  stop("Ignoring tip data names. Number of tip data do not match number of tree tips.")
              }
          }
    }

    ## node data checks
    ## if node.data exist
    if (!all(dim(object@node.data)==0)) {
        ## if we want to use node.names
        if (use.node.names) {
            
            ## check for default names
            if (all(node.names == 1:length(node.names)) 
                || all(node.names == (nTips(object)+1):nEdges(object))) {
                ## no node.names
                if (default.node.names == "fail") {
                    stop("Node data have default names and may not match tree node labels. ",
                         "Consider using the use.node.names=FALSE option.")
                }
                else if (default.node.names == "warn") {
                    warning("Node data have default names and may not match tree node labels. ",
                            "Consider using the use.node.names=FALSE option.")
                }
            }
            
            ## check node names
            ## check for missing or extra node data (relative to tree taxa)
            if (setequal(node.names, object@node.label)) {
                ## names are perfect match - ok
                return(TRUE)
            }
            else {
                ## we know the tree taxa and node.data taxa are not a perfect match
                ## if node.data taxa are subset of tree taxa, check missing.node.data arg and act accordingly
                nodes.in.rownames <- object@node.label %in% node.names
                rownames.in.nodes <- node.names %in% object@node.label 
                missing.data.names <- object@node.label[!nodes.in.rownames]
                missing.data.name.msg <- if (length(missing.data.names)==0) "" else {
                    paste("\n(missing data names: ",
                          paste(missing.data.names,collapse=","),")",sep="")
                }
                extra.data.names <- node.names[!rownames.in.nodes]
                extra.data.name.msg <- if (length(extra.data.names)==0) "" else {
                    paste("\n(extra data names: ",
                          paste(extra.data.names,collapse=","),")",sep="")
                }
                if (!all(nodes.in.rownames)) {
                    ## we know it's not an exact match - we have missing.node.data - take action
                    if (!any(nodes.in.rownames)) {
                        errmsg <- paste("No matches between node data names and tree node labels.",
                                        missing.data.name.msg,extra.data.name.msg)
                        if (missing.node.data == "fail") {
                            stop(errmsg)
                        }
                        else if (missing.node.data == "warn") {
                            warning(errmsg)
                        }
                    }
                    else
                      {
                          errmsg <- paste("Node data names are a subset of tree node labels",
                                          missing.data.name.msg,
                                          extra.data.name.msg)
                          if (missing.node.data == "fail") {
                              stop(errmsg)
                          }

                          else if (missing.node.data == "warn") {
                              warning(errmsg)
                          }
                      }
                    ##else ok
                }
                
                ##if tree taxa are subset of node.data, check extra.node arg and act accordingly
                if (!all(node.names %in% object@node.label)) {
                    ##we know it's not an exact match - we have extra.node.data - take action
                    ##fail
                    errmsg <- paste("Node data names are a superset of tree node labels",
                                    missing.data.name.msg,
                                    extra.data.name.msg)
                    if (extra.node.data == "fail") {
                        stop(errmsg)
                    }
                    ##warn
                    else if (extra.node.data == "warn") {
                        warning(errmsg)
                    }
                    ##else ok
                }
                
                return(TRUE)
            } 
        }
        else
          {
              ##don't use node names or attempt to sort - but check to make sure dimensions match
              if (!(nNodes(object)==dim(object@node.data)[1])) {
                  stop("Ignoring node data names. Number of node data do not match number of tree nodes.")
              }
          }
    }
}

attach_data <- function(object,
                        label.type=c("row.names","column"),
                        label.column=1,
                        use.tip.names=TRUE,
                        use.node.names=FALSE,
                        ...)							 
{
    
    ## assumes data have already been checked by check_data!
    ## name matching default: use row.names of data frame
    label.type = match.arg(label.type)
    if (identical(label.type, "row.names")) {
        tip.names <- row.names(object@tip.data)
        node.names <- row.names(object@node.data)
    }
    else {
        tip.names <- object@tip.data[,label.column]
        node.names <- object@node.data[,label.column]        
    }


    ## for each set of data, take appropriate actions
    
    ## tip data operations:
    ## if tip.data exist
    if (!all(dim(object@tip.data)==0)) {
        ## if we want to use tip.names
        if (use.tip.names) {
            object@tip.data <- object@tip.data[match(object@tip.label,tip.names),,drop=FALSE]
        }
        #tip.names <- object@tip.label
    }
    
    ## node data operations
    if (!all(dim(object@node.data)==0)) {
        ## if we want to use tip.names
        if (use.node.names) {
            object@node.data <- object@node.data[match(object@node.label,node.names),,drop=FALSE]
        }
        #node.names <- object@node.label
    }
    
    return(object)
    
}
