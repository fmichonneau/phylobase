#
# This is the file where to put graphics
# Most functions basically are wrappers for ape
# or ade4 functions
#



##############
# plot phylo4
##############

setGeneric("plot")
setMethod("plot",signature(x="phylo4",y="missing"), function(x,...){
  if(!require(ape)) stop("the ape package is required")
  x <- as(x, "phylo")
  res <- plot(x, ...)
  return(invisible(res))
}) # end plot phylo4



###############
# plot phylo4d
###############
setMethod("plot", signature(x="phylo4d",y="missing"), function(x, type=c("symbols", "squares", "dotchart"), ...){
  if(!require(ade4)) stop("the ade4 package is required")

  if (any(is.na(tdata(x,which="tip")))) {
      warning("dropping NA values from tip data for plot")
      x <- na.omit(x)
  }
  dat <- tdata(x, which="tip")

     
  x <- as(x,"phylog")
  type <- match.arg(type)
  if(ncol(dat)>1 & type=="symbols") type <- "squares" 
  if(ncol(dat)==1 | type=="symbols"){

    res <- symbols.phylog(x, squares=dat[,1], ...)
    
  } else if(type == "squares"){

    res <- table.phylog(dat,x,...)

  } else if(type == "dotchart"){

    res <- dotchart.phylog(x, dat, ...)

  }

  return(invisible(res))
 
})
