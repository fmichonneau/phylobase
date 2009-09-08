################
## subset phylo4
################

setGeneric("subset")
setMethod("subset", "phylo4",
          function(x,tips.include=NULL,tips.exclude=NULL,
                   mrca=NULL,node.subtree=NULL,...) {
              ##  FIXME: could eliminate NULL and make the test
              ## if (!missing) rather than if (!is.null)
            ## (might have to change next line?)
            if (sum(!sapply(list(tips.include,tips.exclude,
                                 mrca,node.subtree),is.null))>1) {
              stop("must specify at most one criterion for subsetting")
            }
            #arglist <- list(...)
            #if (length(arglist)>0) {
            #  warning("unused arguments: ",
            #          paste(names(arglist),collapse=","))
            #}
            kept <- x@tip.label
            dropped <- character(0)
            if (!is.null(tips.include)) {
              if (is.numeric(tips.include)) {
                tips.include <- x@tip.label[tips.include]
              }
              unknown <- setdiff(tips.include,x@tip.label)
              if (length(unknown)>0) {
                warning("unknown tip labels ignored:",
                        paste(unknown,collapse=", "))
                tips.include <- intersect(tips.include,x@tip.label)
              }
              kept <- tips.include
              dropped <- setdiff(x@tip.label,tips.include)
            }
            if (!is.null(tips.exclude)) {
              if (is.numeric(tips.exclude)) {
                tips.exclude <- x@tip.label[tips.exclude]
              }
              unknown <- setdiff(tips.exclude,x@tip.label)
              if (length(unknown)>0) {
                warning("unknown tip labels ignored:",
                        paste(unknown,collapse=", "))
                tips.exclude <- intersect(tips.exclude,x@tip.label)
              }
              dropped <- tips.exclude
              kept <- setdiff(x@tip.label,tips.exclude)
            }
            if (!is.null(node.subtree)) {
              kept <- intersect(x@tip.label,names(descendants(x,node.subtree)))
              dropped <- setdiff(x@tip.label,kept)
            }
            if (!is.null(mrca)) {
              mnode <- MRCA(x,mrca)
              kept <- intersect(x@tip.label,names(descendants(x,mnode)))
              dropped <- setdiff(x@tip.label,kept)
            }
            if (length(kept)<2) {
              stop("0 or 1 tips would remain after subsetting")
            }
            if (length(dropped)==0) return(x)
            return(prune(x, dropped, ...))
          })

###############
# '[' operator
###############
## phylo4
setMethod("[","phylo4",
          function(x, i, j, ..., drop=FALSE) {

              if(missing(i)) i <- TRUE

              oldlab <- tipLabels(x)
              if(is.character(i)){
                  newlab <- i
              } else {
                  newlab <- oldlab[i]
              }
              tip.include <- match(newlab, oldlab)
              res <- subset(x, tips.include=tip.include)

              return(res)
          })




## phylo4d
setMethod("[","phylo4d", 
          function(x, i, j, ..., drop=FALSE) {

              if(missing(i)) i <- TRUE
              if(missing(j)) j <- TRUE

              #### data handling
              ## for now handle only tip data - assumes tip names are good row.names
              tab <- tdata(x, type="tip")[i, j, ..., drop=FALSE]
              
              #### tree handling
              tip.include <- match(row.names(tab), x@tip.label)
              tre <- subset(as(x,"phylo4"), tips.include=tip.include)

              ## result
              res <- phylo4d(x=tre, tip.data=tab)
              
              return(res)
          })

## extract the phylo4 part of phylo4d; relies on implicit coerce method
extractTree <- function(from) {
    as(from, "phylo4")
  }
