require(methods)
require(ape)
         

## accessor functions for all internal bits
## HORRIBLE KLUGE
nTips <- function(x,...)  { }  ## mask ape::nTips
setGeneric("nTips", function(x,...) {
    standardGeneric("nTips")
})
setMethod("nTips","phylo4", function(x,...) {
    ## length(x@tip.label)
    E <- edges(x)
    res <- sum(!E[,2] %in% E[,1])
    return(res)
})
## rm(nTips)

## hack to ensure ape compatibility
setMethod("nTips","ANY", function(x) {
    if (class(x)=="phylo") {
        Ntip(x)
    } else stop(paste("no 'nTips' method available for",
                      deparse(substitute(x)),
                      "(class",class(x),")"))
})

setGeneric("nNodes", function(x) {
    standardGeneric("nNodes")
})
setMethod("nNodes","phylo4", function(x) {
    x@Nnode
})

setGeneric("nEdges", function(x) {
    standardGeneric("nEdges")
})
setMethod("nEdges","phylo4", function(x) {
    nrow(x@edge)
})

setGeneric("edges", function(x,order,...) {
    standardGeneric("edges")
})
setMethod("edges","phylo4", function(x,order,...) {
    x@edge
})

setGeneric("rootEdge", function(x,order,...) {
    standardGeneric("rootEdge")
})
setMethod("rootEdge","phylo4", function(x,order,...) {
    x@root.edge
})

setGeneric("isRooted", function(x) {
    standardGeneric("isRooted")
})

setMethod("isRooted","phylo4", function(x) {
    ## hack to avoid failure on an empty object
    if(nTips(x)==0) return(FALSE)  
    !is.na(x@root.edge) ||  ## root edge explicitly defined
    ## HACK: make sure we find the right "nTips"
    tabulate(edges(x)[, 1])[nTips(x)+1] <= 2
    ## root node (first node after last tip) has <= 2 descendants
    ## FIXME (?): fails with empty tree
})

setGeneric("rootNode", function(x) {
    standardGeneric("rootNode")
})


setMethod("rootNode","phylo4", function(x) {
    if (!isRooted(x)) return(NA)
    if (!is.na(x@root.edge)) stop("FIXME: don't know what to do in this case")
    return(nTips(x)+1)
})

setGeneric("rootNode<-", function(x,value) {
    standardGeneric("rootNode<-")
})

setMethod("rootNode<-","phylo4", function(x,value) {
    stop("not implemented yet")
})

setGeneric("hasEdgeLength", function(x) {
    standardGeneric("hasEdgeLength")
})
setMethod("hasEdgeLength","phylo4", function(x) {
    length(x@edge.length)>0
})

setGeneric("edgeLength", function(x) {
    standardGeneric("edgeLength")
})
setMethod("edgeLength","phylo4", function(x) {
    if (!hasEdgeLength(x)) NULL else x@edge.length
})


setGeneric("hasNodeLabels", function(x) {
    standardGeneric("hasNodeLabels")
})
setMethod("hasNodeLabels","phylo4", function(x) {
    length(x@node.label)>0
})

setGeneric("hasEdgeLabels", function(x) {
    standardGeneric("hasEdgeLabels")
})
setMethod("hasEdgeLabels","phylo4", function(x) {
    length(x@edge.label)>0
})

setGeneric("labels")
setMethod("labels","phylo4", function(object,...) {
    object@tip.label
})
setMethod("labels","phylo4", function(object,which=c("tip","node","allnode"),...) {
    which <- match.arg(which)
    switch(which,tip=object@tip.label,node=object@node.label,
           allnode=c(object@tip.label,object@node.label))
})

setGeneric("labels<-",
           function(object,...,value) {
               standardGeneric("labels<-")
           })
setMethod("labels<-","phylo4", function(object,...,value) {
    if (length(value) != length(object@tip.label))
        stop("Number of tip labels does not match number of tips.")
    object@tip.label <- value
    object
})

setGeneric("nodeLabels", function(x) {
    standardGeneric("nodeLabels")
})
setMethod("nodeLabels","phylo4", function(x) {
    x@node.label
})

setGeneric("nodeLabels<-",
           function(object,...,value) {
               standardGeneric("nodeLabels<-")
           })

setMethod("nodeLabels<-","phylo4", function(object,...,value) {
    object@node.label <- value
    object
})


setGeneric("edgeLabels", function(x) {
    standardGeneric("edgeLabels")
})
setMethod("edgeLabels","phylo4", function(x) {
    x@edge.label
})

setGeneric("edgeLabels<-",
           function(object,...,value) {
               standardGeneric("edgeLabels<-")
           })

setMethod("edgeLabels<-","phylo4", function(object,...,value) {
    object@edge.label <- value
    object
})


## hack to allow access with $
setMethod("$","phylo4",function(x,name) {
    switch(name,
           edge.length=if(!hasEdgeLength(x)) NULL else x@edge.length,
           node.label=if(!hasNodeLabels(x)) NULL else x@node.label,
           root.edge=if(is.na(x@root.edge)) NULL else x@root.edge,
           attr(x,name))
})

## FIXME: implement more checks on this!!
##  do we want to be this permissive?? -- fixed
setMethod("$<-","phylo4",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})

printphylo <- function (x,printlen=6,...) {
    printlen <- max(1,printlen)
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    nb.edge <- length(x$edge.label)
    cat(paste("\nPhylogenetic tree with", nb.tip, "tips and", 
              nb.node, "internal nodes\n"))

    ## print tip labels
    cat("\nTip labels:\n")
    if (nb.tip > printlen) {
        cat(paste("\t", paste(x$tip.label[1:printlen], collapse = ", "), 
                  ", ...\n", sep = ""))
    } else print(x$tip.label)
    
    ## print node labels
    cat("\nNode labels:\n")
    if (nb.node > printlen) {
        cat(paste("\t", paste(x$node.label[1:printlen], collapse = ", "), 
                  ", ...\n", sep = ""))
    } else print(x$node.label)
    
    ## print edge labels
    cat("\nEdge labels:\n")
    if (nb.edge > printlen) {
        cat(paste("\t", paste(x$edge.label[1:printlen], collapse = ", "), 
                  ", ...\n", sep = ""))
    } else print(x$edge.label)

    ## slots
    ##     cat("\nSlots:\n")
    ##     cat(paste("@", names(x)[1:4], sep=""),sep="\t")
    ##     cat("\n")
    ##     cat(paste("@", names(x)[5:7], sep=""),sep="\t")
    ##     cat("\n")
    
    rlab <- if (isRooted(x)) "Rooted"  else "Unrooted"
    cat("\n", rlab, "; ", sep = "")
    blen <- if (hasEdgeLength(x))
      "includes branch lengths"
    else       "no branch lengths"
    cat(blen, "\n\n", sep = "")
}


## hack for print/show 
## from http://tolstoy.newcastle.edu.au/R/e2/devel/06/12/1363.html


##
# Alternative print method for phylo4, showing the contents of the tree data.
##  Not sure if it works for unrooted trees

printphylo4 <- function(x, printall = TRUE){
    if (printall)
      print(as(x, 'data.frame'))
    else print(head(as(x, 'data.frame')))
}

    
setGeneric("print")

#setMethod("print", "phylo4", printphylo)
#setMethod("show", "phylo4", function(object) printphylo(object))
setMethod("print", "phylo4", printphylo4)
setMethod("show", "phylo4", function(object) printphylo4(object))


#################
## summary phylo4
#################
## have to check that x$root.edge is NULL if missing
setMethod("summary","phylo4", function (object, quiet=FALSE) {
    x <- object
    res <- list()
    
    ## build the result object
    res$name <- deparse(substitute(object, sys.frame(-1)))
    res$nb.tips <- length(x$tip.label)
    res$nb.nodes <- x$Nnode
    
    if(!is.null(x$edge.length)){
        res$mean.el <- mean(x$edge.length, na.rm=TRUE)
        res$var.el <- var(x$edge.length, na.rm=TRUE)
        res$sumry.el <- summary(x$edge.length)[-4]
    } else {
        res$mean.el <- NULL
        res$var.el <- NULL
        res$sumry.el <- NULL
    }

    ## polytomies
    if(hasPoly(x)){ # if there are polytomies
        E <- edges(x)
        temp <- tabulate(E[,1])
        degree <- temp[E[,1]] # contains the degree of the ancestor for all edges
        endsAtATip <- !(E[,2] %in% E[,1])
        terminPoly <- (degree>2) & endsAtATip
        internPoly <- (degree>2) & !endsAtATip
        res$degree <- degree
        res$polytomy <- rep("none",nrow(E))
        res$polytomy[terminPoly] <- "terminal"
        res$polytomy[internPoly] <- "internal"
        ## now just keep information about nodes (not all edges)
        nod <- unique(E[,1])
        idx <- match(nod,E[,1])
        res$degree <- res$degree[idx]
        names(res$degree) <- nodeLabels(x)
        res$polytomy <- res$polytomy[idx]
        names(res$polytomy) <- nodeLabels(x)
    }
    
    ## model info
    res$loglik <- attr(x, "loglik")
    res$para <- attr(x, "para")
    res$xi <- attr(x, "xi")
    
    ## if quiet, stop here
    if(quiet) return(invisible(res))
    
    if(!is.null(x$root.edge)){
        cat("  Root edge:", x$root.edge, "\n")
    } else {
        cat("  No root edge.\n")
    }
    ## now, print to screen is !quiet
    cat("\n Phylogenetic tree :", res$name, "\n\n")
    cat(" Number of tips    :", res$nb.tips, "\n")
    cat(" Number of nodes   :", res$nb.nodes, "\n")
    ## cat("  ")
    if(is.null(x$edge.length)) {
        cat(" Branch lengths    : No branch lengths.\n")
    } else {
        cat(" Branch lengths:\n")
        cat("        mean         :", res$mean.el, "\n")
        cat("        variance     :", res$var.el, "\n")
        cat("        distribution :\n")
        print(res$sumry.el)
    }
    if(hasPoly(x)){
        cat("\nDegree of the nodes  :\n")
        print(res$degree)
        cat("\n")
        cat("Types of polytomy:\n")
        print(res$polytomy)
        cat("\n")
    }
    
    if (!is.null(attr(x, "loglik"))) {
        cat("Phylogeny estimated by maximum likelihood.\n")
        cat("  log-likelihood:", attr(x, "loglik"), "\n\n")
        npart <- length(attr(x, "para"))
        for (i in 1:npart) {
            cat("partition ", i, ":\n", sep = "")
            print(attr(x, "para")[[i]])
            if (i == 1)
              next
            else cat("  contrast parameter (xi):", attr(x,"xi")[i - 1], "\n")
        }
    }
    return(invisible(res))
} # end summary phylo4
          ) # end setMethod summary phylo4



setGeneric("tdata", function(x,...) {
    standardGeneric("tdata")
})


setGeneric("hasNodeData", function(x) {
    standardGeneric("hasNodeData")
})

################
## names methods
################
setMethod("names", signature(x = "phylo4"), function(x){
    temp <- rev(names(attributes(x)))[-1]
    return(rev(temp))
})


###################
## Function .genlab
###################
## recursive function to have labels of constant length
## base = a character string
## n = number of labels
.genlab <- function(base,n) {
    if (n<=0) return("")
    s <- seq(length.out=n)
    fw <- max(nchar(as.character(s)))
    numstr <- formatC(s,flag="0",width=fw)
    paste(base,numstr,sep="")
}



## convert from phylo to phylo4
## coerce phylo4d to phylo4 -- on purpose, so no warning
extract.tree <- function(from) {
    phylo4(edge=from@edge,
           edge.length=from@edge.length,
           Nnode=from@Nnode,
           tip.label=from@tip.label)
}


setGeneric("na.omit")
