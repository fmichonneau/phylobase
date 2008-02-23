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
setMethod("plot", signature(x="phylo4d",y="missing"), 
          function(x, treetype=c("phylogram","cladogram"), symbol=c("circles", "squares"), center=TRUE, scale=TRUE, legend=TRUE, grid=TRUE, box=TRUE, show.tip.label=TRUE, show.node.label=TRUE, show.var.label=TRUE, ratio.tree=1/3, font=3, tip.label=x@tip.label, var.label=colnames(x@tip.data), cex.symbol=1, cex.label=1, cex.legend=1, ...){

              
    #### preliminary stuff and checks
    if(!require(ape)) stop("the ape package is required")   
    ## if(ncol(tdata(x,which="tip")) == 0) stop("no data in this phylo4d object")
    
    cex <- par("cex")
    symbol <- match.arg(symbol)
    treetype <- match.arg(treetype)
    
    ## convert the tree into phylo
    tre <- suppressWarnings(as(x,"phylo"))
    tre$node.label <- x@node.label # this should be done by the as(x,"phylo")
    ## plot only tree if no tip data
    if(ncol(tdata(x,which="tip")) == 0) {
         plot.phylo(tre, type=treetype, direction="rightwards", show.tip.label=show.tip.label,
                          show.node.label=show.node.label, cex=cex.label,
                          no.margin=FALSE, x.lim=NULL, y.lim=NULL, ...)
         return(invisible())
     }   

    #### data handling
    ## retrieve data
    dat <- tdata(x, which="tip")
    clas <- lapply(dat,class)
    isNum <- sapply(clas, function(e) e %in% c("integer","numeric"))
    ## keep only numeric data
    dat <- dat[isNum]
    var.label <- var.label[isNum]
   
    ## centring / scaling
    dat <- as.data.frame(scale(dat,center=center,scale=scale))
    
    ## compute bottom margin
    ## ! use inches as units still these won't be changed by plot.phylo
    temp <- var.label[which.max(nchar(var.label))] # longest tip label
    lab.height <- strwidth(temp, units="inches", cex=cex.label) # height required by the longest var label
    lab.height <- lab.height / par("pin")[1] # returned as a fraction of the plot region
    
    #### define plot region
    plotreg <- plotreg0 <- par("plt")
    plotreg.width <- plotreg0[2] - plotreg0[1]
    plotreg.height <- plotreg0[4] - plotreg0[3]
    plotreg[2] <- plotreg[1] + (ratio.tree)*plotreg.width # restrain the width for phylo
    plotreg[3] <- plotreg[3] + plotreg.height*ifelse(show.var.label,lab.height+0.05,0.05) ## add internal vertical margins
    plotreg[4] <- plotreg[4] - plotreg.height*0.05 # add internal vertical margins
    
    #### plot the tree
    par(plt = plotreg)
    plotres <- plot.phylo(tre, type=treetype, direction="rightwards", show.tip.label=FALSE,
                          show.node.label=show.node.label, cex=cex.label,
                          no.margin=FALSE, x.lim=NULL, y.lim=NULL, ...)
    
    #### plot the data
    par(plt=plotreg0)
    cur.usr.width <- par("usr")[2] - par("usr")[1] # beware: par("usr") does not adapt to the new plot region
    usr.width <- cur.usr.width / ratio.tree
    usr.height <- par("usr")[4] - par("usr")[3]
    
    ## x.inset is the space between tree/data and data/tip.labels (in usr units)
    x.inset <- 0.2*cex.symbol * usr.width / par("pin")[1]
    y.inset <- 0.2*cex.symbol * usr.height / par("pin")[2]
    x.base <- plotres$x.lim[2] + x.inset # start plotting from x.base rightwards
    temp <- x@tip.label[which.max(nchar(x@tip.label))] # longest tip label
    lab.width <- strwidth(temp, units="user", cex=cex.label) # compute the width to keep for tip labels
    xrange.data <- c(x.base , (par("usr")[1]+usr.width) - lab.width - 2*x.inset) # plot data within this range

    if(diff(xrange.data) < (x.inset*ncol(dat))) stop("No room left to plot data; please try reducing ratio.tree or cex.label.")

    ## define x and y coordinates
    x.grid <- seq(xrange.data[1],xrange.data[2], length=ncol(dat))
    if(ncol(dat)==1) {x.grid <- mean(c(xrange.data[1],xrange.data[2]))}
    y.grid <- seq(plotres$y.lim[1],plotres$y.lim[2],length=plotres$Ntip)
    temp <- expand.grid(y.grid, x.grid) # here are coordinates for data
    xy.data <- data.frame(x=temp[,2],y=temp[,1])
    
    ## merge data and their coordinates
    alldat <- cbind.data.frame(xy.data, unlist(dat)) 
    fac <- factor(rep(1:ncol(dat), rep(nrow(dat),ncol(dat))))
    alldat <- split(alldat, fac)
    
    ## need to "reboot" the plot region without changing coordinates
    ## seems that box does the job.
    if(box) {box()} else {box(col="transparent")}
    if(grid){
        ## vertical segments
        segments(x0=x.grid, y0=rep(min(y.grid),plotres$Ntip),
                 x1=x.grid, y1=rep(max(y.grid),plotres$Ntip), col="grey")
        ## horizontal segments
        segments(x0=rep(min(x.grid),plotres$Ntip), y0=y.grid,
                 x1=rep(max(x.grid),plotres$Ntip), y1=y.grid, col="grey")
    }

    ## auxiliary function to plot a single variable
    ## max size of a symbol is set to 0.2*cex inches
    ## if changes here, beware to change the 0.15 in x.inset as well
    plotaux <- function(x,y,var,symbol,cex){
        if(any(var[!is.na(var)]<0)) {
            usebw <- TRUE
        } else {
            usebw <- FALSE
        }
        
        if(usebw){
            ispos <- var>0
            fg.col <- rep("black",length(var))
            fg.col[ispos] <- "white"
            bg.col <- rep("white",length(var))
            bg.col[ispos] <- "black"
            
            if(symbol == "squares"){
                symbols(x=x, y=y, squares=abs(var), inches=0.2*cex, fg=fg.col, bg=bg.col, add=TRUE)
            } # end squares
            
            if(symbol == "circles"){
                symbols(x=x, y=y, circles=abs(var), inches=0.2*cex, fg=fg.col, bg=bg.col, add=TRUE)
            } # end circles
            
        } else {

            if(symbol == "squares"){
                symbols(x=x, y=y, squares=var, inches=0.2*cex, fg="white", bg="black", add=TRUE)
            } # end squares
            
            if(symbol == "circles"){
                symbols(x=x, y=y, circles=var, inches=0.2*cex, fg="white", bg="black", add=TRUE)
            } # end circles
        } # end else
        
        if(any(is.na(var))){
            isNA <- is.na(var)
            points(x[isNA],y[isNA],pch=4,cex=cex.symbol)
        }
    } # end plotaux


    ## finally plot the data
    lapply(alldat, function(X) plotaux(X[,1],X[,2],X[,3],symbol,cex.symbol))

    #### plot labels for variables
    if(show.var.label) text(x=x.grid, y=rep(min(y.grid)-1.5*y.inset, ncol(dat)), lab=var.label,
                            adj=1, srt=90, cex=cex.label)

    #### plot tip labels
    x.base <- xrange.data[2] + x.inset
    text(x=rep(x.base,plotres$Ntip), y=1:plotres$Ntip, lab=tip.label, font=font, cex=cex.label, pos=4)


    #### add a legend for symbols
    if(legend){
        leg.var <- alldat[[1]][,3]
        leg.values <- pretty(leg.var,n=4, min.n=1)
        temp <- length(leg.values)
        ## make sure to get maximum 4 symbols
        if(temp>4) {
            leg.values <- leg.values[c(1,2,temp-1,temp)]
        }
        leg.txt <- as.character(leg.values)

        ## temp is a matrix with two columns:
        ## first contains widths of annotations
        ## second contains maximum width of symbols
        temp <- cbind(strwidth(leg.txt,units="user",cex=cex.label*cex.legend) , x.inset*2*cex.legend)
        leg.widths <- apply(temp,1,max)*1.05
        leg.height <- max(strheight(leg.txt, units="user",cex=cex.label*cex.legend))

        ## find basic coordinates
        x.base <- par("usr")[1]+ 0.01*usr.width
        temp <- lab.height * usr.height / (1 - lab.height) ## need to substract temp from par("usr")[3]
        y.base <- par("usr")[3] - temp - y.inset ## to get closer the actual par("usr")[3] !
      
        ## plot annotations
        leg.x <- x.base + leg.widths
        leg.x <- cumsum(leg.x)
        text(leg.x, y.base, leg.txt, cex=cex.label*cex.legend)

        ## plot symbols
        leg.y <- y.base + 2*y.inset*cex.legend
        leg.y <- rep(leg.y,length(leg.x))
        plotaux(leg.x, leg.y, leg.values, symbol, cex.symbol*cex.legend)

        ## FIXME ##
        ## draw a rectangle around the legend
        #rect.size <- c(diff(range(leg.x)) , diff(c(y.base, max(leg.y))) )
        #rect(min(leg.x)- rect.size[1]*0.05,
        #     min(y.base) - rect.size[2]*0.05,
        #     max(leg.x) + rect.size[1]*0.05,
        #     max(y.base) + rect.size[2]*0.05)
    } ## end legend
    
    return(invisible())
}) # end plot phylo4d
