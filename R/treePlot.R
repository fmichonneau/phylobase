
`treePlot` <- function(x, 
                     type = c('phylogram', 'cladogram', 'fan'), 
                     show.tip.label = TRUE,
                     show.node.label = FALSE, 
                     tip.order = NULL,
                     plot.data = is(x, 'phylo4d'),
                     rot = 0,
                     tip.plot.fun = 'bubbles',
                     plot.at.tip = TRUE,
                     edge.color = 'black', 
                     node.color = 'black', # TODO what do with node.color parameter
                     tip.color  = 'black', 
                     edge.width = 1, # TODO line-type modification hack
                     newpage = TRUE,
                     margin = 0.1, # here margin is a precentage of the area
                     ...
            )
{
    ## TODO three dimensional histogram as example, compute values on full dataset
    ## then generate phylo4d object with summary data and plot
    if (!inherits(x, 'phylo4')) stop('treePlot requires a phylo4 or phylo4d object')
    if (!isRooted(x)) stop("treePlot function requires a rooted tree.")
    width  <- height <- (1 - margin)  ## TODO: do these have to be hard-coded?
    type   <- match.arg(type)
    Nedges <- nEdges(x)
    Ntips  <- nTips(x)
    ## TODO remove the false cladogram option?
    if(is.null(edgeLength(x)) || type == 'cladogram') {
        x@edge.length <- rep(1, Nedges)
    }
    xxyy   <- phyloXXYY(x, tip.order)
    x    <- xxyy$x
    tindex <- x@edge[x@edge[, 2] <= Ntips, 2]
    if(type == 'cladogram') {
        xxyy$xx[x@edge[, 2] <= Ntips] <- 1
    }
    
    ## general function for pushing tree subplots
    pushTree <- function(row, col) {
            pushViewport(viewport(layout.pos.row = row, 
                                  layout.pos.col = col,
                                  name = 'treevp'))
                tree.plot(xxyy = xxyy, type = type, 
                    show.tip.label = show.tip.label, 
                    show.node.label = show.node.label, 
                    edge.color = edge.color, node.color = node.color,
                    tip.color = tip.color, edge.width = edge.width, rot = rot)
            upViewport()
    }
    
    # TODO add symbols at the nodes, allow coloirng and sizing downViewport approach?
    # TODO cladogram methods incorrect
    # TODO abstract, make ultrametric? good algorithms for this?
    # TODO for very long plots, alternative margin setting useful
    # call plot.new so that gridBase plots work properly
    # calls to base plot functions need to be cleared w/ par(new = T) which fails
    # if no plot is present TODO perhpas there's a better solution than calling plot.new
    
    ## because we may reoder the tip, we need to update the phy objec
    if(newpage) grid.newpage()
    if(!plot.data) {
        phyplotlayout <- grid.layout(nrow = 1, ncol = 1)
        pushViewport(viewport(width = width, height = height, 
                            layout = phyplotlayout, angle = rot,
                            name = "phyplotlayout"))
            pushTree(row = 1, col = 1)
        upViewport()
        # TODO should return something useful
        return(invisible())
    } else {
        tmin <- min(tdata(x, which = 'tip'), na.rm = TRUE)
        tmax <- max(tdata(x, which = 'tip'), na.rm = TRUE)
        if(!is.function(tip.plot.fun)) {
            if(tip.plot.fun == "bubbles") {
                # use phylobubbles as default
                if(rot != 0) {stop("Rotation of bubble plots not yet implemented")}
                dlabwdth <- max(stringWidth(colnames(x@tip.data))) * 1.2
                phyplotlayout <- grid.layout(nrow = 2, ncol = 2, 
                    heights = unit.c(unit(1, 'null'), dlabwdth), 
                    widths = unit(c(1, 1), c('null', 'null'), list(NULL, NULL)))
                pushViewport(viewport(width = width, height = height, 
                            layout = phyplotlayout, name = 'phyplotlayout'))
                pushViewport(viewport(layout.pos.row = 1:2, layout.pos.col = 2,
                            height = unit(1, 'npc') +
                                            convertUnit(dlabwdth, 'npc'),
                            name = 'bubbleplots', default.units = 'native'))
                    bubout <- phylobubbles(xxyy, ...)
                upViewport()
                pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1,
                            name = 'bubblelegend'))
                    legcir <- seq(bubout$min, bubout$max, length.out = 4)
                    ## print(convertUnit(bubout$bubscale, 'npc', valueOnly = TRUE))
                    ## TODO legend currently does not resize properly
                    legcirS <- legcir * convertUnit(bubout$bubscale, 'inches', valueOnly = TRUE) / bubout$max
                    ccol <- ifelse(legcirS < 0, 'black', 'white')
                    legcirS <- unit(legcirS, 'npc')
                    grid.circle(seq(.2, .8, length.out = length(legcirS)), 0.5, legcirS, gp = gpar(fill = ccol), default.units = 'npc')
                    grid.text(as.character(signif(legcir, digits = 2)), seq(.2, .8, length.out = length(legcir)), 0.1, gp = gpar(cex = 0.75))
                upViewport()
                    pushTree(row = 1, col = 1)
                upViewport()
                return(invisible())
                
            } else if(tip.plot.fun == "density") {
                if(!require(gridBase)) {
                    stop('To plot using base graphics (including the "density"              
                            plot) you need install the "gridBase" package')
                }
                plot.new()
                tmin <- min(tdata(x, which = 'tip'), na.rm = TRUE)
                tmax <- max(tdata(x, which = 'tip'), na.rm = TRUE)
                tip.plot.fun <- function(x,tmin,tmax,...) {
                    # par(omi = c(0,0,0,0))
                    suppressWarnings(par(plt = gridPLT(), new = TRUE))
                    if(!all(is.na(x))) {
                        # hack, set the plotting region to the grid fig region
                        dens <- density(x, na.rm = TRUE)
                        plot.density(dens, xlim = c(tmin, tmax), axes = FALSE,      
                                     mar = c(0,0,0,0), main = "", xlab = "", ylab = "", ...)
                    }
                }
            }
        } else { ## if (is.function(tip.plot.fun))
            ## plot.at.tip <- TRUE
            if (plot.at.tip) {
            datalayout <- grid.layout(ncol = 2,
                                          width = unit(c(1, 1/Ntips), c('null', 'null')) 
                )
            # TODO this is done multiple times, 
            pushViewport(viewport(width = width, height = height, 
                                layout = datalayout, angle = rot,
                                name = 'datalayout'))
            pushViewport(viewport(
                yscale = c(-0.5/Ntips, 1 + 0.5/Ntips), 
                xscale = c(0, 1 + 1/Ntips), 
                layout.pos.col = 1, 
                name = 'data_plots'))
            ## TODO should plots float at tips, or only along edge?
            hc <- convertY(unit(1/Ntips, 'snpc'), 'npc')
            for(i in 1:Ntips) {
                pushViewport(viewport(
                    y = xxyy$yy[x@edge[, 2] == i],
                    x = 1 + 1/(2 * Ntips), # xxyy$xx[phy@edge[, 2] == i], 
                    height = hc, 
                    width = hc, 
                    # default.units = 'native', 
                    name = paste('data_plot', i),
                    just = "center",
                    angle = -rot
                    ))
                    #grid.rect()
                    vals = t(tdata(x, which = 'tip')[i, ])
                    if (!all(is.na(vals))) tip.plot.fun(vals,tmin,tmax,...)
                upViewport()
            }
            pushTree(row = 1, col = 1)
            upViewport()
            upViewport()
        } else { ## plot by data column
            ## !plot.at.tip
            ## stop("not implemented yet")
            nvars <- length(x@tip.data) ## FIXME
            datalayout <- grid.layout(ncol = 2)  ## set default widths equal for this kind of plot
            ## width = unit(c(1, 1/Ntips), c('null', 'null'))
            # TODO this is done multiple times, 
            pushViewport(viewport(width = width, height = height, 
                                layout = datalayout, angle = rot,
                                name = 'datalayout'))
            pushViewport(viewport(
                yscale = c(-0.5/Ntips, 1 + 0.5/Ntips), 
                xscale = c(0, 1 + 1/Ntips), 
                layout.pos.col = 2, 
                name = 'data_plots'))
            ## TODO should plots float at tips, or only along edge?
            hc <- convertY(unit(1/nvars, 'snpc'), 'npc')
            for(i in 1:nvars) {
                pushViewport(viewport(
                    x = i/nvars,   ## xxyy$yy[phy@edge[, 2] == i],
                    y =  0.5, ## 1 + 1/(2 * Ntips), # xxyy$xx[phy@edge[, 2] == i], 
                    height = 1, ## hc, 
                    width = 1/nvars, ## hc, 
                    # default.units = 'native', 
                    name = paste('data_plot', i),
                    just = "center",
                    angle = -rot
                    ))
                    #grid.rect()
                    vals = tdata(x)[,i]
                    if (!all(is.na(vals))) tip.plot.fun(vals, tmin, tmax, ...)
                upViewport()
            }
            upViewport()
            pushTree(row = 1, col = 1)
            upViewport()
        }
            return(invisible())
        } ## if (is.function(tip.plot.fun))
    }
}

tree.plot <- function(xxyy, type, show.tip.label, show.node.label, edge.color, 
                        node.color, tip.color, edge.width, rot) 
{
    # TODO switch to phylobase abstractions
    phy <- xxyy$phy
    Nedges   <- nEdges(phy)
    Ntips    <- nTips(phy)
    tindex <- phy@edge[phy@edge[, 2] <= Ntips, 2]
    eindex <- match(phy@edge[,2], xxyy$phy.orig@edge[,2])
    segs <- segs(XXYY = xxyy)

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
    nindex <- sort(eindex[phy@edge[, 2] > Ntips], index.return = TRUE)$ix
    if(length(node.color) != length(nindex)) {
        node.color <- rep(node.color, length.out = length(nindex))
    }
    node.color <- node.color[nindex]

    if(show.tip.label) {
        ## calculate several lab dimesisions
        ## labw    -- a vector of string widths
        ## adjlabw -- the max width for adjusting the size of viewports
        ## laboff  -- a vector of half string widths for 
        ## offsetting center justified labels, handy for vp rotation 
        labw    <- stringWidth(phy@tip.label)
        adjlabw <- max(labw) + unit(0.02, 'npc')
        ## print(foo <<- adjlabw)
        laboff  <- labw * 0.5
        treelayout <- grid.layout(nrow = 1, ncol = 2,
            widths = unit.c(unit(1, 'null', NULL), convertUnit(adjlabw, 'inches'))
            )
        tindex <- phy@edge[phy@edge[, 2] <= Ntips, 2]
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
        layout = treelayout, name = 'treelayout'))
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
    upViewport()
    if(show.tip.label) {
        pushViewport(viewport(layout.pos.col = 2,
            name = 'tiplabelvp'))
        labtext <- grid.text(
            phy@tip.label[tindex], 
            x = convertUnit(laboff[tindex], 'native', valueOnly = TRUE) + 0.02,
            y = xxyy$yy[phy@edge[, 2] %in% tindex], rot = lrot,
            default.units = 'native', name = 'tiplabels',
            just = 'center', gp = gpar(col = tip.color[tindex])
        )
        upViewport()
    }
    # TODO probably want to be able to adjust the location of these guys
    if(show.node.label) {
        pushViewport(viewport(layout = treelayout, layout.pos.col = 1))
            rty <- mean(xxyy$yy[phy@edge[, 1] == Ntips + 1], name = 'nodelabelvp')
        labtext <- grid.text(
            phy@node.label, 
            x = c(0, xxyy$xx[phy@edge[, 2] > Ntips][nindex]), 
            y = c(rty, xxyy$yy[phy@edge[, 2] > Ntips][nindex]), 
            default.units = 'npc', name = 'nodelabels', rot = -rot,
            just = 'center', gp = gpar(col = node.color[nindex])
        )
        upViewport()
    }
    upViewport()
    # grobTree(vseg, hseg, labtext)
}


phyloXXYY <- function(phy, tip.order = NULL) {
    ## initalize the output
    Nedges <- nEdges(phy)
    phy.orig <- phy
    Ntips  <- nTips(phy)
    xxyy = list(
        yy = rep(NA, Nedges), 
        xx = numeric(Nedges), 
        ## record the order that nodes are visited in
        traverse = NULL) 
    if(is.null(edgeLength(phy))) {
        # TODO there should be an abstraction for assigning branch lengths
        stop('Phylogeny has no branch lengths, cannot calculate x coordinates')
    }
    
    # TODO tip ordering should be dealt with at a higher level
    # if(!is.null(tip.order)) { 
    #     yy[which(phy@edge[, 2] == tip.order)] <- seq(
        ## TODO perhaps we want to use match here?
        ## 0, 1, length.out = Ntips) 
    # } else {
        ## reoder the phylo and assign even y spacing to the tips
        phy <- reorder(phy, 'pruningwise')
        xxyy$yy[phy@edge[, 2] <= Ntips] <- seq(
            0, 1, length.out = Ntips
        )
    # }
    
    ## a recurvise preorder traversal 
    ## node  -- initalized to be root, is the starting point for the traversal
    ## phy   -- the phylogeny
    ## xxyy  -- the list initalized below that holds the output
    ## prevx -- the sum of ancestral branch lengths
    calc.node.xy <- function(node, phy, xxyy, prevx = 0) {
        ## if node == root node, and there is no root edge set get descendants
        ## and set index to NULL index is used for indexing output
        if(any(phy@edge[, 2] == node) == FALSE) {
            decdex <- which(phy@edge[, 1] == node)
            index <- NULL
            ## if root node start at x = 0
            newx <- 0
        } else {
            ## non-root node behavior
            ## get row in edge matrix corresponding to node, get descendants
            index <- which(phy@edge[, 2] == node)
            decdex <- which(phy@edge[, 1] == phy@edge[index, 2])
            ## non-root node x location 
            newx <- xxyy$xx[index] <- phy@edge.length[index] + prevx
            ## if the x value is already set we are at a tip and we return
            if(!is.na(xxyy$yy[index])) { return(xxyy) }
        }
        for(i in phy@edge[decdex, 2]) {
            ## for each decendant call the function again
            xxyy <- calc.node.xy(i, phy, xxyy, newx)
        }
        if(!is.null(index)) {
            ## set y value by averaging the decendants
            xxyy$yy[index] <- mean(xxyy$yy[decdex])
        }
        ## TODO performance improvement here? rely on above ordering?
        ## keep track of the nodes and order we visited them
        xxyy$traverse <- c(xxyy$traverse, phy@edge[decdex, 2]) 
        xxyy
    }
    ## call function for the first time
    xxyy <- calc.node.xy(Ntips + 1, phy, xxyy)
    ## scale the x values
    xxyy$xx <- xxyy$xx / max(xxyy$xx)
    # TODO return an index vector instead of a second phy object
    c(xxyy, phy = list(phy), phy.orig = list(phy.orig))
}

segs <- function(XXYY) {
    phy <- XXYY$phy
    treelen <- rep(NA, nEdges(phy) + 1)
    segs <- list(v0x = treelen, v0y = treelen, v1x = treelen, v1y = treelen,
                 h0x = treelen, h0y = treelen, h1x = treelen, h1y = treelen)
    troot <- nTips(phy) + 1

    get.coor <- function(node, segs) {
        if(any(phy@edge[, 2] == node) == FALSE) {
            #root
            decdex <- which(phy@edge[, 1] == node)
            index <- length(treelen)
            segs$v0y[index] <- segs$v1y[index] <- NA
            segs$v0x[index] <- segs$v1x[index] <- NA
            segs$h0y[index] <- segs$h1y[index] <- NA
            segs$h0x[index] <- segs$h1x[index] <- NA
            segs$v0x[decdex] <- segs$v1x[decdex] <- segs$h0x[decdex] <- 0            
            segs$v0y[decdex] <- mean(XXYY$yy[decdex])            
        } else {
            #not root
            index <- which(phy@edge[, 2] == node)
            segs$h1x[index] <- XXYY$xx[index]
            segs$h0y[index] <- segs$h1y[index] <- segs$v1y[index] <- XXYY$yy[index]
            if(!any(phy@edge[, 1] == node)) {
                return(segs)
            }
            decdex <- which(phy@edge[, 1] == phy@edge[index, 2])
            segs$v0x[decdex] <- segs$v1x[decdex] <- segs$h0x[decdex] <- XXYY$xx[index]            
            segs$v0y[decdex] <- mean(XXYY$yy[decdex])           
        }
        
        for(i in phy@edge[decdex, 2]) {
            segs <- get.coor(i, segs)
        }
        segs
    }
    get.coor(troot, segs)
}

phylobubbles <- function(XXYY, square = FALSE, grid = TRUE) {
    ## TODO add legend command
    ## tys   -- tip y coordinates
    ## nVars -- number of traits/characters
    ## maxr  -- maximum circle radius, based on nVars or nTips
    
    phy <- XXYY$phy
    
    # tip y coordinates
    tys <- XXYY$yy[phy@edge[, 2] <= nTips(phy)]
    
    tipdata <- tdata(phy, which = "tip")
    nVars <- ncol(tipdata) # number of bubble columns
    
    maxr <- ifelse(ncol(tipdata) > nTips(phy), 1 / ncol(tipdata), 1 / nTips(phy))
    tipdata <- apply(tipdata, 2, 
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
    dnas <- is.na(tipdata)
    naxs <- naxs[dnas]
    nays <- nays[dnas]
    ## set the NA points to zero so that grid.circle doesn't crash
    tipdata[is.na(tipdata)] <- 0
    
    ## get label widths
    tiplabwidth  <- max(stringWidth(phy@tip.label))
    datalabwidth <- max(stringWidth(colnames(tipdata)))
    
    ## 2x2 layout -- room at the bottom for data labels, and legend
    bublayout <- grid.layout(nrow = 2, ncol = 2,
        widths  = unit.c(unit(1, 'null', NULL), tiplabwidth), 
        heights = unit.c(unit(1, 'null', NULL), datalabwidth * 1.2))
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
        sqedge <- unit(unlist(tipdata), 'snpc')
        grid.rect(x = xrep, y = yrep, 
            width = sqedge, 
            height = sqedge, 
            gp=gpar(fill = ccol))
    } else {
        ## plot bubbles
        grid.circle(xrep, yrep, r = unlist(tipdata), gp = gpar(fill = ccol))
    }
    # catch a value for scaling other output
    bubscale <- convertUnit(unit(max(tipdata, na.rm = TRUE), 'npc'), 'inches')
    upViewport()
    
    ## push view ports for tip and data labels fixed locations
    pushViewport(viewport( 
        name = 'bubble_tip_labels', 
        layout = bublayout, 
        layout.pos.col = 2, 
        layout.pos.row = 1
    ))
    grid.text(phy@tip.label, 0.1, tys, just = 'left')
    upViewport()
    pushViewport(viewport( 
        name = 'bubble_data_labels', 
        layout = bublayout, 
        layout.pos.col = 1, 
        layout.pos.row = 2
    ))
    grid.text(colnames(tipdata), xpos, .65, rot = 90, just = 'right')

    upViewport(2)
    # to make a nice legend, return the biggest smallest and a scaling factor
    # translate the scale of the current vp to a fixed value
    ## ensure the min is not a zero (or NA) that's replaced by a zero
    ## print(convertUnit(bubscale, 'inches', valueOnly = TRUE))
    return(list(max = max(tipdata, na.rm = TRUE), 
                min = min(tipdata[tipdata != 0], na.rm = TRUE),
                has.na = length(naxs) > 0,
                bubscale = bubscale))
}

# setGeneric("treePlot", useAsDefault = treePlot)
# setMethod("treePlot", signature = c('phylo4', 'phylo4d'), treePlot)

# gridbasefun <- function(f, naked = TRUE, scale = TRUE) {
#     function(x, tmin, tmax, ...) {
#         require(gridBase)
#         op <- par()
#         if (naked) {
#             par(ann = FALSE, mar = rep(0, 4))
#             ## this could take a bit of hacking
#             ## to work nicely in general -- too bad
#             ## par(ann=FALSE) doesn't work in general
#             ## Could set main="", xlab="", ylab="", axes=FALSE
#             ## but this will break 
#         }
#         ## this must be the *last* par() call
#         suppressWarnings(par(plt = gridPLT(), new = TRUE)) 
#         if(!all(is.na(x))) {
#             f(x, xlim = c(tmin, tmax), ...)
#         }
#     }
# }
# 
# setGeneric('plot', useAsDefault = treePlot)
setMethod('treePlot', signature = c('phylo4', 'phylo4d'), treeplot)
