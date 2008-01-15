################
## plot phylo4d
################
#setMethod("plot", signature(x="phylo4d",y="missing"), 
plottemp <- function(x, symbol=c("circles", "squares"), usebw=TRUE, center=TRUE, scale=TRUE, legend=TRUE, grid=TRUE, box=TRUE, show.tip.label=TRUE, show.node.label=TRUE, show.var.label=TRUE, ratio.tree=1/3, font=3, tip.label=x@tip.label, var.label=colnames(x@tip.data), cex.symbol=1, cex.label=par("cex"), cex.legend=1, ...){

    #### preliminary stuff and checks
    if(ncol(tdata(x,which="tip")) == 0) stop("no data in this phylo4d object")
    if (any(is.na(tdata(x,which="tip")))) {
        warning("dropping NA values from tip data for plot")
        x <- na.omit(x)
    }

    cex <- par("cex")
    
    symbol <- match.arg(symbol)

    ## convert the tree into phylo
    tre <- suppressWarnings(as(x,"phylo"))
    tre$node.label <- x@node.label # this should be done by the as(x,"phylo")

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
    
   
    #### define plot region
    plotreg <- plotreg0 <- par("plt")
    plotreg.width <- plotreg0[2] - plotreg0[1]
    plotreg.height <- plotreg0[4] - plotreg0[3]
    plotreg[2] <- plotreg[1] + (ratio.tree)*plotreg.width ## restrain the width for phylo
    plotreg[3] <- plotreg[3] + plotreg.height*ifelse(show.var.label,0.20,0.05) ## add internal vertical margins
    plotreg[4] <- plotreg[4] - plotreg.height*0.05 ## add internal vertical margins
    
    #### plot the tree
    par(plt = plotreg)
    plotres <- plot(tre, direction="rightwards", show.tip.label=FALSE,
                    show.node.label=show.node.label, cex=cex.label, ...)
    
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
    plotaux <- function(x,y,var,symbol,usebw,cex){
        if(any(var<0)) usebw <- TRUE
        
        if(usebw){
            ispos <- var>0
            fg.col <- rep("white",length(var))
            fg.col[ispos] <- "black"
            bg.col <- rep("black",length(var))
            bg.col[ispos] <- "white"
            
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
    } # end plotaux


    ## plot the data
    lapply(alldat, function(X) plotaux(X[,1],X[,2],X[,3],symbol,usebw,cex.symbol))

    ## add legend for variables
    if(show.var.label) text(x=x.grid, y=rep(min(y.grid)-1.5*y.inset, ncol(dat)), lab=var.label,
                            adj=1, srt=90, cex=cex.label)

    #### plot tip labels
    x.base <- xrange.data[2] + x.inset
    text(x=rep(x.base,plotres$Ntip), y=1:plotres$Ntip, lab=tip.label, font=font, cex=cex.label, pos=4)

    return(invisible())
} # end plot phylo4d





#####
## tests
##
#example("phylo4d")

obj1 <- obj2 <- obj3 <- phylo4d(as(tree.owls,"phylo4"),data.frame(wing=1:4,color=factor(c("b","w","b","b")), tail=runif(4)*10), use.tip.names=FALSE)

obj2@tip.data <- as.data.frame(obj2@tip.data[,1])
obj3@tip.data <- cbind(obj1@tip.data,obj2@tip.data)

plottemp(obj1)
plottemp(obj2,box=FALSE)
plottemp(obj3)
par(mar=rep(.5,4))
plottemp(obj3,box=FALSE,cex.sym=1.2,cex.la=.8)


library(ade4)
data(mjrochet)
temp <- as(read.tree(text=mjrochet$tre),"phylo4")
obj <- phylo4d(x=temp,tip.data=mjrochet$tab)
obj
plot(obj)
