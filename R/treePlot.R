`treePlot` <- function(phy, 
                     type = c('phylogram', 'cladogram', 'fan'), 
                     show.tip.label = TRUE,
                     show.node.label = FALSE, 
                     tip.order = NULL,
                     plot.data = is(phy, 'phylo4d'),
                     rot = 0,
                     tip.plot.fun = 'bubbles',
                     plot.at.tip = TRUE,
                     edge.color = 'black', 
                     node.color = 'black', # TODO what do with node.color parameter
                     tip.color  = 'black', 
                     edge.width = 1, # TODO line-type modification hack
                     newpage = TRUE,
                     margins = c(1.1, 1.1, 1.1, 1.1), # number of lines, same as par(mar)
                     ...
            )
{
    ## TODO three dimensional histogram as example, compute values on full dataset
    ## then generate phylo4d object with summary data and plot

    ## TODO factors not handled in data plots
    ## TODO add symbols at the nodes, allow coloirng and sizing downViewport approach?
    ## TODO cladogram methods incorrect
    ## because we may reoder the tip, we need to update the phy objec

    if (!inherits(phy, 'phylo4')) stop('treePlot requires a phylo4 or phylo4d object')
    if (!isRooted(phy)) stop("treePlot function requires a rooted tree.")
    if (plot.data && !hasTipData(phy)) {
        warning("tree has no tip data to plot")
        plot.data <- FALSE
    }
    
    if(newpage) grid.newpage()
    type   <- match.arg(type)
    Nedges <- nEdges(phy)
    Ntips  <- nTips(phy)
    
    if(!is.null(tip.order)) {
        if(length(tip.order) != Ntips) {stop('tip.order must be the same length as nTips(phy)')}
        if(is.numeric(tip.order)) {
            tip.order <- tip.order
        } else {
            if(is.character(tip.order)) {
                tip.order <- match(tip.order, tipLabels(phy))
            }
        }
    }
    
    ## TODO remove the false cladogram option?
    if(!hasEdgeLength(phy) || type == 'cladogram') {
        edgeLength(phy) <- rep(1, Nedges)
    }
    xxyy   <- phyloXXYY(phy, tip.order)
    phy    <- xxyy$phy
    pedges <- edges(phy)
    tindex <- pedges[pedges[, 2] <= Ntips, 2][tip.order]
    if(type == 'cladogram') {
        xxyy$xx[pedges[, 2] <= Ntips] <- 1
    }
    
    ## plotViewport is a convience function that provides margins in lines
    pushViewport(plotViewport(margins=margins))
    
    if(!plot.data) {
        plotOneTree(xxyy, type, show.tip.label, show.node.label, edge.color, 
                                node.color, tip.color, edge.width, rot)
    } else {
        if(!is.function(tip.plot.fun)) {
            if(tip.plot.fun == "bubbles") {
                phylobubbles(
                    type = type, 
                    show.node.label = show.node.label, 
                    tip.order = tip.order, 
                    rot = 0, 
                    edge.color = edge.color, 
                    node.color = node.color, # TODO what do with node.color parameter
                    tip.color  = tip.color, 
                    edge.width = edge.width, # TODO line-type modification hack
                    newpage = TRUE, 
                    ..., XXYY = xxyy 
                )
            } else {
                stop(paste(tip.plot.fun, 'is neither a function or a recognized plot type'))
            }
        } else { ## from -- if(tip.plot.fun == "bubbles")
            ## plot.at.tip <- TRUE
            if (plot.at.tip) {
                tip.data.plot(
                    xxyy = xxyy, 
                    type = type, 
                    show.tip.label = show.tip.label, 
                    show.node.label = show.node.label, 
                    tip.order = tip.order, 
                    rot = 0, 
                    tip.plot.fun = tip.plot.fun, 
                    edge.color = edge.color, 
                    node.color = node.color, # TODO what do with node.color parameter
                    tip.color  = tip.color, 
                    edge.width = edge.width, # TODO line-type modification hack
                    newpage = TRUE, 
                    ... 
                )
                return(invisible())
            } ## if (plot.at.tip)
        } ## else
    } ## else
    upViewport() # margins
}

plotOneTree <- function(xxyy, type, show.tip.label, show.node.label, edge.color, 
                        node.color, tip.color, edge.width, rot) 
{
    # TODO switch to phylobase abstractions
    phy    <- xxyy$phy
    Nedges <- nEdges(phy)
    Nnodes <- nNodes(phy)
    Ntips  <- nTips(phy)
    pedges <- edges(phy)
    tindex <- pedges[pedges[, 2] <= Ntips, 2]
    eindex <- match(pedges[,2], edges(xxyy$phy.orig)[,2])
    segs   <- xxyy$segs

    ## TODO check that colors are valid?
    if(length(edge.color) != Nedges) {
        edge.color <- rep(edge.color, length.out = Nedges)
    }
    edge.color <- edge.color[eindex]
    
    if(length(edge.width) != Nedges) {
        edge.width <- rep(edge.width, length.out = Nedges)
    }
    edge.width <- edge.width[eindex]

    ## TODO check that colors are valid?
    if(length(node.color) != Nnodes) {
        node.color <- rep(node.color, length.out = Nnodes)
    }

    if(show.tip.label) {
        ## calculate several lab dimesisions
        ## labw    -- a vector of string widths
        ## adjlabw -- the max width for adjusting the size of viewports
        ## laboff  -- a vector of half string widths for 
        ## offsetting center justified labels, handy for vp rotation 
        labw    <- stringWidth(tipLabels(phy))
        adjlabw <- max(labw) + unit(0.1, 'inches')
        laboff  <- labw * 0.5 + unit(0.1, 'inches')
        ## print(foo <<- laboff)
        treelayout <- grid.layout(nrow = 1, ncol = 2,
            widths = unit.c(unit(1, 'null', NULL), convertUnit(adjlabw, 'inches'))
            )
        tindex <- pedges[pedges[, 2] <= Ntips, 2]
        if(length(tip.color) != Ntips) {
            tip.color <- rep(tip.color, length.out = Ntips)
        }
        # keep labels horizontal unless plot is upwards or downwards
        lrot <- ifelse(rot %% 360 %in% c(90, 270), 0, -rot)
    } else {
        treelayout <- grid.layout(nrow = 1, ncol = 1)
    }
    # grid.show.layout(treelayout)
    pushViewport(viewport(
        x = 0.5, y = 0.5, 
        width = 1, height = 1, 
        layout = treelayout, angle = rot, name = 'treelayout'))
    pushViewport(viewport(
        layout.pos.col = 1, 
        name = 'tree'))
    if (type == "fan") {
        dseg <- grid.segments( # draws diag lines
            x0 = segs$v0x, y0 = segs$v0y, 
            x1 = segs$h1x, y1 = segs$h1y, 
            name = "diag", gp = gpar(col = edge.color, lwd = edge.width))     
    } else {
        vseg <- grid.segments( # draws vertical lines
            x0 = segs$v0x, y0 = segs$v0y, 
            x1 = segs$v1x, y1 = segs$v1y, 
            name = "vert", gp = gpar(col = edge.color, lwd = edge.width)) 
        hseg <- grid.segments( # draws horizontal lines
            x0 = segs$h0x, y0 = segs$h0y, 
            x1 = segs$h1x, y1 = segs$h1y, 
            name = "horz", gp = gpar(col = edge.color, lwd = edge.width))
    }
    upViewport() # tree
    if(show.tip.label) {
        pushViewport(viewport(layout.pos.col = 1,
            name = 'tiplabelvp'))
        labtext <- grid.text(
            tipLabels(phy)[tindex], 
            x = unit(xxyy$xx[pedges[, 2] %in% tindex], "native") + laboff[tindex],
            y = xxyy$yy[pedges[, 2] %in% tindex], rot = lrot,
            default.units = 'native', name = 'tiplabels',
            just = 'center', gp = gpar(col = tip.color[tindex])
        )
        upViewport() #tiplabelvp
    }
    # TODO probably want to be able to adjust the location of these guys
    if(show.node.label) {
        pushViewport(viewport(layout = treelayout, layout.pos.col = 1, name = 'nodelabelvp'))
            theLabels <- nodeLabels(phy)
            # don't plot NAs
            theLabels[is.na(theLabels)] <- ""
        labtext <- grid.text(
            theLabels, 
            x = c(xxyy$xx[pedges[, 2] > Ntips]), 
            y = c(xxyy$yy[pedges[, 2] > Ntips]), 
            default.units = 'npc', name = 'nodelabels', rot = -rot,
            just = 'center', gp = gpar(col = node.color)
        )
        upViewport() #nodelabelvp
    }
    upViewport() # treelayout
    # grobTree(vseg, hseg, labtext)
}

phyloXXYY <- function(phy, tip.order = NULL) 
{
    phy.orig <- phy
    ## initalize the output
    phy    <- reorder(phy, 'preorder')
    pedges <- edges(phy)
    Nedges <- nrow(pedges) ## TODO switch to the accessor once stablized
    Ntips  <- nTips(phy)
    tips <- pedges[, 2] <= Ntips
    if(!is.null(tip.order)) {
        tip.order <- match(tip.order, pedges[, 2][tips])
    }
    xx <- numeric(Nedges)
    yy <- numeric(Nedges)

    treelen <- rep(NA, nEdges(phy))
    segs <- list(v0x = treelen, v0y = treelen, v1x = treelen, v1y = treelen,
                 h0x = treelen, h0y = treelen, h1x = treelen, h1y = treelen)
    
    ## Set root x value to zero and calculate x positions
    xx[1] <- 0
    segs$v0x[1] <- segs$v1x[1] <- segs$h0x[1] <- 0 
    edge1   <- as.integer(pedges[,1])
    edge2   <- as.integer(pedges[,2])
    edgeLen <- edgeLength(phy)
    edgeLen[is.na(edgeLen)] <- 0
    edgeLen <- as.numeric(edgeLen)
    nedges  <- as.integer(nEdges(phy))
    segsv0x <- as.numeric(rep.int(0, Nedges))
    xPos <- .C("phyloxx", edge1, edge2,
            edgeLen, nedges, xx, segsv0x)
    ## browser()
    xx <- xPos[[5]]
    segs$v0x <- xPos[[6]]
    ## test1 <- function() {
    ##     for (i in edge[, 2]) {
    ##         dex <- edge[, 1] == i
    ##         cur <- edge[, 2] == i
    ##         xx[dex] <- phy@edge.length[dex] + xx[cur]
    ##         segs$v0x[dex] <- xx[cur]
    ##     }
    ##     return(list(segs=segs, xx=xx))
    ## }
    ## test1out <- test1()
    ## segs <- test1out$segs
    ## xx   <- test1out$xx

    ## Set y positions for terminal nodes and calculate remaining y positions
    if(!is.null(tip.order)) {
        yy[tips][tip.order] <- seq(0, 1, length = Ntips)
    } else {
        yy[tips] <- seq(0, 1, length = Ntips)
    }
    segs$h0y[tips] <- segs$h1y[tips] <- yy[tips]
    segs$v1y[tips] <- segs$v0y[tips] <- yy[tips]
    placeHolder <- function() {
        for(i in rev((Ntips + 1):nEdges(phy))) {
            dex <- pedges[, 1] == i
            cur <- pedges[, 2] == i
            yy[cur] <- segs$v0y[dex] <- mean(yy[dex])
        }
        return(list(segs=segs, yy=yy))
    }
    placeHolder2 <- function() {
        for(i in rev((Ntips + 1):nEdges(phy))) {
            cur <- pedges[, 2] == i
            dex <- pedges[, 1] == i
            yy[cur] <- segs$v0y[dex] <- mean(yy[dex])
        }
        return(list(segs=segs, yy=yy))
    }

    yPos <- placeHolder()
    segs <- yPos$segs
    yy   <- yPos$yy

    ## edgeLen[is.na(edgeLen)] <- 0
    ## edgeLen <- as.numeric(edgeLen)
    ## ntips   <- as.integer(nTips(phy))
    ## yy      <- as.numeric(yy)
    ## segsv0y <- as.numeric(yy)
    ## browser()
    ## yPos <- .C("phyloyy", edge1, edge2,
    ##         ntips, nedges, yy, segsv0y)

    segs$h0y <- segs$h1y <- segs$v1y <- yy

    ## scale the x values so they range from 0 to 1
    Xmax <- max(xx)
    segs$v0x <- segs$v0x / Xmax
    xx <- xx / Xmax
    
    segs$h1x <- xx
    segs$v1x <- segs$h0x <- segs$v0x 
    
    # TODO return an index vector instead of a second phy object
    list(xx = xx, yy = yy, phy = phy, phy.orig = phy.orig, segs = segs)
}


.bubLegendGrob <- function(tipdata, tipdataS) {
    grob(tipdata=tipdata, tipdataS=tipdataS, cl='bubLegend')
}

drawDetails.bubLegend <- function(x, ...) {
    ## number of bubbles in legend
    leglen  <- 4
    ## the raw tip data
    tipdata <- x$tipdata
    ## the tip data as scaled for bubble plot
    ts <- x$tipdataS
    ## return to the bubble plot viewport to get properly scaled values
    ## this relies on having well named unique viewports
    seekViewport("bubble_plots")
        ## retreive the min and max non-zero bubbles as numerics not units
        bubrange <- convertUnit(
                    unit(c(min(ts[ts != 0], na.rm=TRUE), max(ts[ts != 0], na.rm=TRUE)), "native"), 
                    "mm", valueOnly=TRUE)
    seekViewport("bubblelegend")
    ## grid.rect()
    ## Generate the sequence of legend bubble sizes and convert to grid mm units
    legcirS  <- unit(seq(bubrange[1], bubrange[2], length.out=leglen), "mm")
    ## get the corresponding sequence of actual data values
    legcir   <- seq(min(tipdata[tipdata != 0], na.rm=TRUE), 
                    max(tipdata[tipdata != 0], na.rm=TRUE), length.out=leglen)
    ccol     <- ifelse(legcir < 0, 'black', 'white')

    leftedge <- abs(convertUnit(legcirS[1], 'npc', valueOnly=TRUE)) + 0.1
    xloc     <- seq(leftedge, 0.5, length.out=leglen)
    textsp   <- convertUnit(max(abs(legcirS)), axisFrom="y", axisTo="y", 'npc', valueOnly=TRUE)
    strsp    <- convertUnit(unit(1, "strheight", "TTT"), axisFrom="y", 'npc', valueOnly=TRUE)
    grid.circle(x=xloc, y=0.9 - textsp - strsp, r=legcirS, gp = gpar(fill=ccol), default.units = 'npc')
    grid.text(as.character(signif(legcir, digits = 2)), 
                x=xloc, y=0.75 - 2 * textsp - strsp, 
                gp=gpar(cex=0.75), 
                default.units='npc'
    )
}

phylobubbles <- function(type = type,
                        place.tip.label = "right", 
                        show.node.label = show.node.label, 
                        tip.order = NULL,
                        rot = 0,
                        edge.color = edge.color, 
                        node.color = node.color, # TODO what do with node.color parameter
                        tip.color  = tip.color, 
                        edge.width = edge.width, # TODO line-type modification hack
                        newpage = TRUE,
                        ..., 
                        XXYY, square = FALSE, grid = TRUE)
{
    ## TODO add legend command
    ## tys   -- tip y coordinates
    ## nVars -- number of traits/characters
    ## maxr  -- maximum circle radius, based on nVars or nTips
    if(rot != 0) {stop("Rotation of bubble plots not yet implemented")}
    lab.right <- ifelse(place.tip.label %in% c("right", "both"), TRUE, FALSE)
    lab.left  <- ifelse(place.tip.label %in% c("left", "both"), TRUE, FALSE)


    phy     <- XXYY$phy
    tmin    <- min(tdata(phy, type = 'tip'), na.rm = TRUE)
    tmax    <- max(tdata(phy, type = 'tip'), na.rm = TRUE)
    tipdata <- tdata(phy, type = "tip")[nodeId(phy,"tip"),,drop=FALSE]
    nVars   <- ncol(tipdata) # number of bubble columns
    pedges  <- edges(phy)

    dlabwdth <- max(stringWidth(colnames(tipdata))) * 1.2
    if(convertWidth(dlabwdth, 'cm', valueOnly=TRUE) < 2) {dlabwdth <- unit(2, 'cm')}
    phyplotlayout <- grid.layout(nrow = 2, ncol = 2, 
        heights = unit.c(unit(1, 'null'), dlabwdth), 
        widths = unit(c(1, 1), c('null', 'null'), list(NULL, NULL)))
    pushViewport(viewport(layout = phyplotlayout, name = 'phyplotlayout'))
    pushViewport(viewport(layout.pos.row = 1:2, layout.pos.col = 2,
                height = unit(1, 'npc') +
                                convertUnit(dlabwdth, 'npc'),
                name = 'bubbleplots', default.units = 'native'))

    # tip y coordinates
    tys <- XXYY$yy[pedges[, 2] <= nTips(phy)]
    
    maxr <- ifelse(ncol(tipdata) > nTips(phy), 1 / ncol(tipdata), 1 / nTips(phy))
    tipdataS <- apply(tipdata, 2, 
                    function(x) (maxr * x) / max(abs(x), na.rm = TRUE))
    if(nVars == 1) {
        xpos <- 0.5
    } else {
        xpos <- seq(0 + maxr + 0.02, 1 - maxr - 0.02, length.out = nVars)
    }

    ## rep coordinates for filling a matrix columnwise
    xrep <- rep(xpos, each = length(tys))
    yrep <- rep(tys, nVars)
    ## color bubbles 
    ccol <- ifelse(tipdata < 0, 'black', 'white')
    
    ## generate matrices of every x and y by filling the repd value columnwise
    ## then subset for datapoints that are NA
    naxs <- matrix(xrep, ncol = nVars)
    nays <- matrix(yrep, ncol = nVars)
    dnas <- is.na(tipdataS)
    naxs <- naxs[dnas]
    nays <- nays[dnas]
    ## set the NA points to zero so that grid.circle doesn't crash
    tipdataS[is.na(tipdataS)] <- 0
    
    ## get label widths
    if(lab.right) {
        tiplabwidth  <- max(stringWidth(tipLabels(phy)))
    } else {tiplabwidth <- unit(0, 'null', NULL)}

    ## 2x2 layout -- room at the bottom for data labels, and legend
    bublayout <- grid.layout(nrow = 2, ncol = 2,
        widths  = unit.c(unit(1, 'null', NULL), tiplabwidth), 
        heights = unit.c(unit(1, 'null', NULL), dlabwdth))
    pushViewport(viewport(
        x = 0.5, y = 0.5, 
        width = 0.95, height = 1, 
        layout = bublayout, name = 'bublayout'
    ))
    pushViewport(viewport( 
        name = 'bubble_plots', 
        layout = bublayout, 
        layout.pos.col = 1, 
        layout.pos.row = 1
    ))
    if(grid) {
        ## draw light grey grid behind bubbles
        grid.segments(x0 = 0,   x1 = 1, 
                      y0 = tys, y1 = tys, gp = gpar(col = 'grey'))
        grid.segments(x0 = xpos, x1 = xpos, 
                      y0 = 0,    y1 = 1, gp = gpar(col = 'grey'))
    }    
    if (length(naxs) > 0) {
        ## if ther are missing values plot Xs
        grid.points(naxs, nays, pch = 4)
    }
    
    if(square) {
        ## alternative to circles
        ## to keep the squares square, yet resize nicely use the square npc
        sqedge <- unit(unlist(tipdataS), 'snpc')
        grid.rect(x = xrep, y = yrep, 
            width = sqedge, 
            height = sqedge, 
            gp=gpar(fill = ccol))
    } else {
        ## plot bubbles
        grid.circle(xrep, yrep, r = unlist(tipdataS), gp = gpar(fill = ccol))
    }
    upViewport()
    
    ## push view ports for tip and data labels fixed locations
    if(lab.right) {
        pushViewport(viewport( 
            name = 'bubble_tip_labels', 
            layout = bublayout, 
            layout.pos.col = 2, 
            layout.pos.row = 1
        ))
        tt <- tipLabels(phy) # phy@tip.label 
        grid.text(tt, 0.1, tys, just = 'left')
        upViewport()
    }
    pushViewport(viewport( 
        name = 'bubble_data_labels', 
        layout = bublayout, 
        layout.pos.col = 1, 
        layout.pos.row = 2
    ))
    ## ideas, for nicer sizing of the data labels
    ## data.label.space <- convertX(unit(1, 'npc'), "points", valueOnly = TRUE)
    ## data.label.fontsize <- data.label.space / ncol(tipdata)
    ## , gp=gpar(fontsize=data.label.fontsize))
    ## offset the data labels from the bottom bubble
    datalaboffset <- convertUnit(unit(15, "mm"), 'npc', valueOnly=TRUE)
    grid.text(colnames(tipdata), xpos, 1-datalaboffset, rot = 90, just = 'right')

    upViewport(3)
    pushViewport(viewport(layout.pos.row=2, layout.pos.col=1,
                name='bubblelegend'))
    yyy <- .bubLegendGrob(tipdata, tipdataS)
    grid.draw(yyy)
    upViewport()
    
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1,
                name = 'tree'))
        plotOneTree(XXYY, type, show.tip.label=lab.left, show.node.label, edge.color, 
                                node.color, tip.color, edge.width, rot)
    upViewport(2)
    
    # to make a nice legend, return the biggest smallest and a scaling factor
    # translate the scale of the current vp to a fixed value
    ## ensure the min is not a zero (or NA) that's replaced by a zero
    ## print(convertUnit(bubscale, 'inches', valueOnly = TRUE))
    ## return(list(max = max(tipdata, na.rm = TRUE), 
    ##             min = min(tipdata[tipdata != 0], na.rm = TRUE),
    ##             has.na = length(naxs) > 0,
    ##             bubscale = bubscale))
}

tip.data.plot <- function(
                     xxyy, 
                     type = c('phylogram', 'cladogram', 'fan'), 
                     show.tip.label = TRUE,
                     show.node.label = FALSE, 
                     tip.order = NULL,
                     rot = 0, 
                     tip.plot.fun = grid.points, 
                     edge.color = 'black', 
                     node.color = 'black', # TODO what do with node.color parameter
                     tip.color  = 'black', 
                     edge.width = 1, # TODO line-type modification hack
                     ...)    
{
    phy    <- xxyy$phy
    pedges <- edges(phy)
    Ntips  <- nTips(phy)
    datalayout <- grid.layout(ncol = 2, width = unit(c(1, 1/Ntips), c('null', 'null')))
    # TODO this is done multiple times, 
    pushViewport(viewport(layout = datalayout, angle = rot,
                        name = 'datalayout'))
    pushViewport(viewport(
        yscale = c(-0.5 / Ntips, 1 + 0.5 / Ntips), 
        xscale = c(0, 1 + 1 / Ntips), 
        layout.pos.col = 1, 
        name = 'data_plots'))
    ## TODO should plots float at tips, or only along edge?
    hc <- convertY(unit(1 / Ntips, 'snpc'), 'npc')
    for(i in 1:Ntips) {
        pushViewport(viewport(
            y = xxyy$yy[pedges[, 2] == i],
            x = 1 + 1 / (2 * Ntips), # xxyy$xx[phy@edge[, 2] == i], 
            height = hc, 
            width = hc, 
            # default.units = 'native', 
            name = paste('data_plot', i),
            just = "center",
            angle = -rot
            ))
            #grid.rect()
            tvals <- tdata(phy, type = 'tip')[nodeId(phy,'tip'), , drop=FALSE]
            vals = t(tvals[i, ])
            if (!all(is.na(vals))) tip.plot.fun(vals, ...)
        upViewport() # loop viewports
    }
    plotOneTree(xxyy, type, show.tip.label, show.node.label, edge.color, 
                            node.color, tip.color, edge.width, rot)    
    upViewport(2) ## data_plot & datalayout
}

# phyloStripchart <- function()

setGeneric('plot')
setMethod('plot', signature(x='phylo4', y='missing'), function(x, y, ...) {
    treePlot(x, ...)
})

