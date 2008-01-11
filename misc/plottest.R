################
## plot phylo4d
################
#setMethod("plot", signature(x="phylo4d",y="missing"), function(x, symbol=c("circles", "squares"), center=TRUE,

plottemp <- function(x, symbol=c("circles", "squares"), center=TRUE, scale=TRUE, legend=TRUE, grid=TRUE, show.tip.label=TRUE, show.node.label=FALSE, ratio.tree=1/3, cex.symbol=1, cex.label=par("cex"), cex.legend=1, ...){
   
    if (any(is.na(tdata(x,which="tip")))) {
        warning("dropping NA values from tip data for plot")
        x <- na.omit(x)
    }

    ## retrieve data
    dat <- tdata(x, which="tip")

    ## keep only numeric data

    #
    ## TO DO 
    #

    ## convert data to a list
    dat <- as.list(dat)
    
    cex <- list(...)$cex
    if(is.null(cex)) cex <- par("cex")
     
    ## tre <- suppressWarnings(as(x,"phylo"))
    tre <- as(x,"phylo")
    
    symbol <- match.arg(symbol)

    ## define plot region
    plotreg <- plotreg0 <- par("plt")
    plotreg.width <- plotreg0[2] - plotreg0[1]
    plotreg[2] <- plotreg[1] + (ratio.tree + 0.01)*plotreg.width
    
    ## plot the tree
    par(plt = plotreg)
    plotres <- plot(tre, direction="rightwards", show.tip.label=FALSE)
    
    ## plot the data
    par(plt=plotreg0)
    usr.width <- par("usr")[2] - par("usr")[1]
    x.base <- plotres$x.lim[2] + usr.width*0.05 # start plotting from x.base rightwards

    tip.labels <- x@tip.label
    temp <- tip.labels[which.max(nchar(tip.labels))] # longest tip label
    lab.width <- strwidth(temp, units="user", cex=cex.label) # compute the width to keep for tip labels

    xrange.data <- c(x.base , par("usr")[2] - lab.width - usr.width*0.05) # plot data within this range

    x.data <- seq(xrange.data[1],xrange.data[2], length=ncol(dat))
    y.data <- seq(plotres$y.lim[1],plotres$y.lim[2],length=plotres$Ntip)
    xy.data <- expand.grid(x.data, y.data) # here are coordinates for data

    if(grid){
        abline(v=x.data,h=y.data, col="grey")
    }
    
    if(symbol == "squares"){
        symbols(x=xy.data[1], y=xy.data[2], squares=as.numeric(dat), inches=0.5, add=TRUE)
    }

    if(symbol == "circles"){
        symbols(x=xy.data[1], y=xy.data[2], circles=as.numeric(dat), inches=0.5, add=TRUE)
           
    }

    ## plot tip labels

    
    return(invisible())
}) # end plot phylo4d
