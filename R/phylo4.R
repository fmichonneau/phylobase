require(methods)
require(ape)

setOldClass("phylo")
## setOldClass("multi.tree") ## obsolete
setOldClass("multiPhylo")

setClass("phylo4",
         representation(edge="matrix",
                        edge.length="numeric",
                        Nnode="integer",
                        node.label="character",
                        tip.label="character",
                        edge.label="character",
                        root.edge="numeric"),
         prototype=list(edge=matrix(nrow=0,ncol=2,dimname=list(NULL,c("ancestor","descendent"))),
           edge.length=numeric(0),
           Nnode=as.integer(0),
           tip.label=character(0),
           node.label=character(0),
           edge.label=character(0),
           ## check?
           ##           node.label = as.character(1:Nnode),
           root.edge=as.numeric(NA)),
         validity=check_phylo4)
         
###################################
## phylo4d class
## extend: phylo with data
setClass("phylo4d",
         representation(tip.data="data.frame",
                        node.data="data.frame"),
         ##                        edgedata="data.frame"),
         prototype = list( tip.data = data.frame(NULL),
           node.data = data.frame(NULL) ),
         ##all.data = data.frame(NULL) ),
         validity = function(object) {
             ## FIXME: finish this by intercepting FALSE, char string, etc.
             check1 <- check_data(object)
             check2 <- check_phylo4(object)          
         },                   
         contains="phylo4")
         

## accessor functions for all internal bits
## HORRIBLE KLUGE
nTips <- function(x,...)  { }  ## mask ape::nTips
setGeneric("nTips", function(x,...) {
    standardGeneric("nTips")
})
setMethod("nTips","phylo4", function(x,...) {
    length(x@tip.label)
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

setGeneric("RootEdge", function(x,order,...) {
    standardGeneric("RootEdge")
})
setMethod("RootEdge","phylo4", function(x,order,...) {
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

setGeneric("EdgeLength", function(x) {
    standardGeneric("EdgeLength")
})
setMethod("EdgeLength","phylo4", function(x) {
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

setGeneric("labels<-",
           function(object,...,value) {
               standardGeneric("labels<-")
           })

setMethod("labels<-","phylo4", function(object,...,value) {
    object@tip.label <- value
    object
})


setGeneric("NodeLabels", function(x) {
    standardGeneric("NodeLabels")
})
setMethod("NodeLabels","phylo4", function(x) {
    x@node.label
})

setGeneric("NodeLabels<-",
           function(object,...,value) {
               standardGeneric("NodeLabels<-")
           })

setMethod("NodeLabels<-","phylo4", function(object,...,value) {
    object@node.label <- value
    object
})


setGeneric("EdgeLabels", function(x) {
    standardGeneric("EdgeLabels")
})
setMethod("EdgeLabels","phylo4", function(x) {
    x@edge.label
})

setGeneric("EdgeLabels<-",
           function(object,...,value) {
               standardGeneric("EdgeLabels<-")
           })

setMethod("EdgeLabels<-","phylo4", function(object,...,value) {
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
setAs(from='phylo4',to='data.frame',
      def = function(from) {
      	x <- from
	ancestor <- x@edge[,1]
	node <- x@edge[,2]
	root <- unique(ancestor[!ancestor %in% node])
	int.node <- c(root, unique(ancestor[ancestor %in% node]))
        tip <- node[!(node %in% ancestor)]
	n.tip <- length(tip)
        n.int <- length(int.node)
        node <- c(root, node)
        if (length(ancestor)>0) ancestor <- c(NA, ancestor)
        branch.length <- c(x@root.edge, x@edge.length)
        if (length(branch.length) == 1) branch.length <- rep("", n.tip+n.int)
        if (print.species <- !(is.null(x@node.label) & is.null(x@tip.label)))
          { 
              nl <- x@node.label       
              if (is.null(nl)) nl <-  rep("", n.int)   # phylo4 has a node.label for the root?
              tl <- x@tip.label
              if (is.null(tl)) tl <-  rep("", n.tip)
              species.name <- c(nl, tl)
          } else species.name <- NULL
        if (length(root)==0) {
            node.type <- c(rep("internal", n.int), rep("tip", n.tip))
        }  else node.type <- c("root", rep("internal", n.int-1), rep("tip", n.tip))
        
        return(data.frame(species.name, node, ancestor, branch.length, node.type))
    })

printphylo4 <- function(x, printall = TRUE){
    if (printall)
      print(as(x, 'data.frame'))
    else print(head(as(x, 'data.frame')))
}

setAs(from='phylo4d', to='data.frame',
      function(from) {

    as(from, "phylo4") -> tree  	 # get tree
    as(tree, "data.frame") -> t_df   # convert to data.frame
    tdata(from, "allnode") -> dat               # get data
    t_df$order <- rownames(t_df)     # save roworder of tree

    merge(t_df, dat, by.x="species.name", by.y="row.names", all.x=TRUE, sort=FALSE) -> tdat
                                     # merged tree, data, but mixed up

    order(tdat$order) -> o           # ordering to get back original roworder                              
    tdat[o,] -> tdat                 # original order of tree
    rownames(tdat) <- tdat$order     # fix scrambled row numbers
    return(tdat[,!colnames(tdat) %in% "order"]) # drop "order"
})
    
setGeneric("print")

#setMethod("print", "phylo4", printphylo)
#setMethod("show", "phylo4", function(object) printphylo(object))
setMethod("print", "phylo4", printphylo4)
setMethod("show", "phylo4", function(object) printphylo4(object))
setMethod("print", "phylo4d", printphylo4)
setMethod("show", "phylo4d", function(object) printphylo4(object))


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
        res$polytomy <- rep("no poly.",nrow(E))
        res$polytomy[terminPoly] <- "terminal poly."
        res$polytomy[internPoly] <- "internal poly."
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
        cat(" Branch lengths      : No branch lengths.\n")
    } else {
        cat(" Branch lengths:\n")
        cat("        mean         :", res$mean.el, "\n")
        cat("        variance     :", res$var.el, "\n")
        cat("        distribution :\n")
        print(res$sumry.el)
    }
    if(hasPoly(x)){
        cat("\nDegree of the nodes  :", res$degree, "\n")
        cat("Polytomies at the nodes:", res$polytomy, "\n")
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
setMethod("tdata","phylo4d", function(x,which=c("tip","node","allnode"),...) {
    which <- match.arg(which)
    if (which=="allnode") {
        namesmatch <- all(colnames(x@tip.data)==colnames(x@node.data))
        classmatch <- all(sapply(x@tip.data,class)==sapply(x@node.data,class))
        if (!(classmatch && namesmatch)) stop("Node and tip columns do not match, access tip and node data separately")
    }
    switch(which,tip=x@tip.data,node=x@node.data,
           allnode=rbind(x@tip.data,x@node.data))
    ##         edge=x@edgedata)
})



## setMethod("summary", "phylo4d", function(object){
##   x <- object
##   tdata(x, "tip") -> tips
##   tdata(x, "allnode") -> allnodes
##   cat("Phylogenetic tree with", nTips(x), " species and", nNodes(x), "internal nodes\n\n")
##   cat("  Tree plus data object of type:", class(x), "\n")
##   cat("  Species Names                :", labels(x), "\n")
##   if (hasEdgeLength(x)){ 
##     cat("  Has Branch Lengths (first 10):", EdgeLength(x)[1:min(length(EdgeLength(x)),10)], "\n")
##   } 
##   cat("  Rooted                       :", isRooted(x), "\n\n\n")
##  
##   cat("\nComparative data\n")
##   if (nrow(tips)>0) 
##     {
##       cat("\nTips: data.frame with", nTips(x), "species and", ncol(tips), "variables \n")
##       print(summary(tips))
##     }
##   if (nrow(allnodes)>0) 
##     {
##       cat("\nNodes: data.frame with", nEdges(x), "species and internal nodes and", ncol(allnodes), "variables \n")                  ## May have to fix once  Node=Edge issue is settled
##       print(summary(allnodes))
##     }
##   
## }) # end summary phylo4d
## 

## Alternative phylo4d summary method, using phylo4 summary
## Marguerite Butler & Peter Cowan
setMethod("summary", "phylo4d", function(object){
    x <- object

    summary(as(object, "phylo4"))

    tdata(object, "tip") -> tips
    tdata(object, "node") -> nodes

    cat("\nComparative data:\n")
    if (nrow(tips) > 0) 
      {
          cat("\nTips: data.frame with", nTips(object), "taxa and", ncol(tips), "variables \n\n")
          print(summary(tips))
      }else {cat('\nObject contains no tip data.')}

    if (nrow(nodes) > 0) 
      {
          cat("\nNodes: data.frame with", nNodes(object), "internal nodes and", ncol(nodes), "variables \n\n")                  ## May have to fix once  Node=Edge issue is settled
          print(summary(nodes))
      } else {cat('\nObject contains no node data.\n')}

}) # end summary phylo4d

## extend: phylo with model fit (???)
## hacked with logLik attribute from ape, but otherwise not done




setClass("multiPhylo4",
         representation(phylolist="list",
                        tree.names="character"),
         prototype = list(phylolist=list(),
           tree.names=character(0)))

setClass("multiPhylo4d",
         representation(tip.data="data.frame"),
         contains="multiPhylo4")

################
## show phylo4d
################
##
setMethod("show", "phylo4d", function(object){
    x <- object

    cat("\n##Comparative data##\n")
    ##  print tree
    cat("\n#Tree#\n")
    printphylo(x)

    ## print traits
    cat("\n#Traits#\n")
    cat("\ntip.data: data.frame containing", ncol(tdata(x,"tip")), "traits for", nrow(tdata(x,"tip")),"tips" )
    cat("\nnode.data: data.frame containing", ncol(tdata(x,"node")), "traits for", nrow(tdata(x,"node")),"nodes" )

    cat("\n")
}) # end summary phylo4d

## ?? setMethod("print", "phylo4", o)


setGeneric("hasNodeData", function(x) {
    standardGeneric("hasNodeData")
})
setMethod("hasNodeData","phylo4d", function(x) {
    nrow(x@node.data)>0
})

setMethod("NodeLabels<-","phylo4d", function(object,...,value) {
    object@node.label <- value
    rownames(object@node.data) <- value
    object
})

setMethod("labels<-","phylo4d", function(object,...,value) {
    object@tip.label <- value
    rownames(object@tip.data) <- value
    object
})


################
## names methods
################
setMethod("names", signature(x = "phylo4"), function(x){
    temp <- rev(names(attributes(x)))[-1]
    return(rev(temp))
})

setMethod("names", signature(x = "phylo4d"), function(x){
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



#####################
## phylo4 constructor
#####################
##
## TEST ME . wait for validity check
##
phylo4 <- function(edge, edge.length=NULL, tip.label=NULL, node.label=NULL,
                   edge.label=NULL, root.edge=NULL,...){
    ## edge
    mode(edge) <- "integer"
    if(any(is.na(edge))) stop("NA are not allowed in edge matrix")
    if(ncol(edge)>2) warning("the edge matrix has more than two columns")
    edge <- as.matrix(edge[,1:2])
    colnames(edge) <- c("ancestor","desendent")
    
    ## edge.length
    if(!is.null(edge.length)) {
        if(!is.numeric(edge.length)) stop("edge.length is not numeric")
        edge.length <- edge.length
    } else {
        edge.length <- as.numeric(NULL)
    }

    ## tip.label
    ntips <- sum(tabulate(edge[,1]) == 0)
    if(is.null(tip.label)) {
        tip.label <- .genlab("T",ntips)
    } else {
        if(length(tip.label) != ntips) stop("the tip labels are not consistent with the number of tips")
        tip.label <- as.character(tip.label)
    } 

    ## node.label
    nnodes <- sum(tabulate(edge[,1]) > 0)
    if(is.null(node.label)) {
        node.label <- .genlab("N",nnodes)
    } else {
        if(length(node.label) != nnodes) stop("the node labels are not consistent with the number of nodes")
    } 

    ## edge.label
    ## an edge is named by the descendant
    if(is.null(edge.label)) {
        edge.label <- paste("E", edge[,2], sep="")
    } else {
        if(length(edge.label) != nrow(edge)) stop("the edge labels are not consistent with the number of edges")
    } 

    ## root.edge - if no root edge lenth provided, set to a numeric NA
    if(is.null(root.edge)) root.edge <- as.numeric(NA)
    ##if(!is.null(root.edge)) {
    ##    if(!round(root.edge)==root.edge) stop("root.edge must be an integer")
    ##    root.edge <- as.integer(root.edge)
    ##    if(root.edge > nrow(edge)) stop("indicated root.edge do not exist")
    ##} else {
    ##    root.edge <- as.integer(NA)
    ##}
    
    ## fill in the result
    res <- new("phylo4")
    res@edge <- edge
    res@edge.length <- edge.length
    res@Nnode <- nnodes
    res@tip.label <- tip.label
    res@node.label <- node.label
    res@edge.label <- edge.label
    res@root.edge <- root.edge

    ## check_phylo4 will return a character string if object is
    ##  bad, otherwise TRUE
    if (is.character(checkval <- check_phylo4(res))) stop(checkval)
    return(res)
}




######################
## phylo4d constructor
######################
## TEST ME 
## '...' recognized args for data are tipdata and nodedata.
## other recognized options are those known by the phylo4 constructor
##

## generic
setGeneric("phylo4d", function(x, ...) { standardGeneric("phylo4d")} )

## first arg is a phylo4
setMethod("phylo4d", c("phylo4"), function(x, tip.data=NULL, node.data=NULL, all.data=NULL, ...){

    if(!check_phylo4(x)) stop("invalid phylo4 object provided in x")
    
    res <- new("phylo4d")
    res@edge <- x@edge
    res@edge.length <- x@edge.length
    res@Nnode <- x@Nnode
    res@tip.label <- x@tip.label
    res@node.label <- x@node.label
    res@edge.label <- x@edge.label
    res@root.edge <- x@root.edge

### handle a which argument
    which.dat <- match.arg(list(...)$"which", c("tip","node","all"))

    ## handle data
    if(all(is.null(c(tip.data, node.data, all.data)))) {
        stop("no data provided; please use phylo4 class")
    }

    ## convert vector to data.frames
    if(is.vector(tip.data)) tip.data <- as.data.frame(tip.data)
    if(is.vector(node.data)) node.data <- as.data.frame(node.data)
    if(is.vector(all.data)) all.data <- as.data.frame(all.data)
    
    if(!is.null(all.data)){
        if(!is.data.frame(all.data)) stop("all.data must be a data.frame")
        tip.data <- all.data[1:nTips(x) , , drop=FALSE]
        node.data <- all.data[-(1:nTips(x)) , , drop=FALSE]
    }

    ## now at least one data.frame is provided
    if(is.null(tip.data)) tip.data <- data.frame(NULL)
    if(is.null(node.data)) node.data <- data.frame(NULL)
    if(!is.data.frame(tip.data)) stop("tip.data must be a data.frame")
    if(!is.data.frame(node.data)) stop("node.data must be a data.frame")
    
    res@tip.data <- tip.data
    res@node.data <- node.data

    check_data(res, ...)
    res <- attach_data(res,...)
    return(res)  
})

## first arg is a matrix of edges
setMethod("phylo4d", c("matrix"), function(x, tip.data=NULL, node.data=NULL, all.data=NULL, ...){
    tree <- phylo4(edge=x,...)
    res <- phylo4d(tree, tip.data, node.data, all.data, ...)
    return(res)
})

## first arg is a phylo
setMethod("phylo4d", c("phylo"), function(x, tip.data=NULL, node.data=NULL, all.data=NULL, ...){
    tree <- as(x, "phylo4")
    res <- phylo4d(tree, tip.data, node.data, all.data, ...)
    return(res)
})

## convert from phylo to phylo4
setAs("phylo","phylo4",
      function(from,to) {
          newobj <- phylo4(from$edge, from$edge.length,
                           from$tip.label,
                           node.label=from$node.label,
                           edge.label=from$edge.label, ## ???
                           root.edge=from$root.edge)
          attribs = attributes(from)
          attribs$names <- NULL
          knownattr <- c("logLik","order","origin","para","xi")
          known <- names(attribs)[names(attribs) %in% knownattr]
          unknown <- names(attribs)[!names(attribs) %in% c(knownattr,"class","names")]
          if (length(unknown)>0) {
              warning(paste("unknown attributes ignored: ",unknown,collapse=" "))
          }
          for (i in known) attr(newobj,i) <- attr(from,i)
          newobj
      })

setAs("phylo","phylo4d",
      function(from,to) {
          phylo4d(as(from,"phylo4"),tip.data=data.frame())
      })

setAs("multiPhylo4","multiPhylo",
      function(from,to) {
          newobj <- new("multiPhylo4",
                        phylolist=lapply(from,as,to="phylo4"))
      })

setAs("multiPhylo4d","multiPhylo",
      function(from,to) {
          newobj <- new("multiPhylo4d",
                        phylolist=lapply(from,as,to="phylo4"),
                        tree.names=names(from),
                        tip.data=data.frame())
      })

setAs("multiPhylo","multiPhylo4",
      function(from,to) {
          y <- lapply(as,from@phylolist,to="phylo")
          names(y) <- from@tree.names
          if (nrow(from@tip.data)>0) warning("discarded tip data")
          class(y) <- "multiPhylo"
          y
      })

setAs("phylo4","phylo",
      function(from,to) {
          y <- list(edge=from@edge,
                    edge.length=from@edge.length,
                    Nnode=from@Nnode,
                    tip.label=from@tip.label,
                    node.label=from@node.label)
          class(y) <- "phylo"
          if(length(y$edge.length) == 0) y$edge.length <- NULL
          if(length(y$node.label) == 0) y$node.label <- NULL
          if (!is.na(from@root.edge)) y$root.edge <- from@root.edge
          y
      })

## coerce phylo4d to phylo4 -- on purpose, so no warning
extract.tree <- function(from) {
    phylo4(edge=from@edge,
           edge.length=from@edge.length,
           Nnode=from@Nnode,
           tip.label=from@tip.label)
}

setAs("phylo4d","phylo",
      function(from,to) {
          y <- list(edge=from@edge,
                    edge.length=from@edge.length,
                    Nnode=from@Nnode,
                    tip.label=from@tip.label)
          class(y) <- "phylo"
          if(length(y$edge.length) == 0) y$edge.length <- NULL
          if(length(y$node.label) == 0) y$node.label <- NULL
          if (!is.na(from@root.edge)) y$root.edge <- from@root.edge
         
          warning("losing data while coercing phylo4d to phylo")
          y
      })



####################
## as(phylo4,phylog)
####################
setOldClass("phylog")
setAs("phylo4","phylog", function(from, to){
    if(!require(ade4)) stop("the ade4 package is required")
    x <- as(from,"phylo")
    x <- write.tree(x,file="")
    x <- newick2phylog(x, add.tools=FALSE)
    return(x)
})

## FIXME: doesn't deal with missing node data
##   (don't even know how that should be done in this case)
setGeneric("na.omit")
setMethod("na.omit", "phylo4d",
          function(object, ...) {
            tipdata <- tdata(object,"tip")
            na.names <- rownames(tipdata)[!complete.cases(tipdata)]
            prune(object,tip=na.names)
          })
