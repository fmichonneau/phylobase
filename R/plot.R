#
# This is the file where to put graphics
# Most functions basically are wrappers for ape
# or ade4 functions
#



###############
## plot phylo4
###############

setGeneric("plot")
setMethod("plot",signature(x="phylo4",y="missing"), function(x,...){
    if(!require(ape)) stop("the ape package is required")
    x <- as(x, "phylo")
    plot(x, ...)
    return(invisible())
}) # end plot phylo4



################
## plot phylo4d
################
setMethod("plot", signature(x="phylo4d",y="missing"), function(x, type=c("symbols", "squares", "dotchart"),
                                          show.tip.label=TRUE, show.node.label=FALSE, cex.label=par("cex"),
                                          cex.symbol=1, ...){
    if(!require(ade4)) stop("the ade4 package is required")
    
    if (any(is.na(tdata(x,which="tip")))) {
        warning("dropping NA values from tip data for plot")
        x <- na.omit(x)
    }
    if (show.node.label) warning("DANGER: node labels may be reordered from original tree")
    dat <- tdata(x, which="tip")

    cex <- list(...)$cex
    if(is.null(cex)) cex <- par("cex")
     
    tre <- suppressWarnings(as(x,"phylog"))
    type <- match.arg(type)
    if(ncol(dat)>1 & type=="symbols") type <- "squares" 
    if(ncol(dat)==1 | type=="symbols"){
        
        symbols.phylog(tre, squares=dat[,1], csize=cex.symbol, ...)
        
    } else if(type == "squares"){
        arglist <- list(df=dat, phylog=tre, labels.nod=rev(x$node.label),
                        clabel.row=as.numeric(show.tip.label)*cex.label,
                        clabel.nod=as.numeric(show.node.label)*cex.label,
                        csize=max(2,cex.symbol),
                        ...)

        do.call(table.phylog, arglist)
        
    } else if(type == "dotchart"){
        arglist <- list(phylog=tre, values=dat, labels.nod=rev(x$node.label),
                        clabel.row=as.numeric(show.tip.label)*cex.label,
                        clabel.nod=as.numeric(show.node.label)*cex.label,
                        ...)
           
        do.call(dotchart.phylog, arglist)
        
    }
    
    return(invisible())
}) # end plot phylo4d
