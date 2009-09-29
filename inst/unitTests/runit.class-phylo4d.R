#
# --- Test class-phylo4d.R ---
#

# create ape::phylo version of a simple tree for testing
nwk <- "((t1:0.1,t2:0.2)n7:0.7,(t3:0.3,(t4:0.4,t5:0.5)n9:0.9)n8:0.8)n6:0.6;"
tr <- read.tree(text=nwk)

# create analogous phylo4 object with a full complement of valid slots
ancestor <- as.integer(c(6,7,7,6,8,0,8,9,9))
descendant <- as.integer(c(7,1,2,8,3,6,9,4,5))
edge <- cbind(ancestor, descendant)
nid.tip <- 1:5
nid.int <- 6:9
nid.all <- c(nid.tip, nid.int)
lab.tip <- paste("t", nid.tip, sep="")
lab.int <- paste("n", nid.int, sep="")
lab.all <- c(lab.tip, lab.int)
elen <- descendant/10
elab <- paste("e", ancestor, descendant, sep="-")
phy <- phylo4(x=edge, tip.label=lab.tip, node.label=lab.int,
    edge.length=elen, edge.label=elab)

# create altered version such that each slot is out of order with
# respect to all others; methods should be able to handle this
phy.alt <- phy
phy.alt@tip.label <- rev(phy@tip.label)
phy.alt@node.label <- rev(phy@node.label)
phy.alt@edge <- phy@edge[c(6:9, 1:5), ]
phy.alt@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy.alt@edge.label <- phy@edge.label[c(8:9, 1:7)]

# create data to add to phylo4 to create phylo4d, but with data rows out
# of order
set.seed(1)
nid.tip.r <- sample(nid.tip)
nid.int.r <- sample(nid.int)
nid.all.r <- sample(c(nid.tip, nid.int))
allDt <- data.frame(a=letters[nid.all.r], b=10*nid.all.r)
tipDt <- data.frame(c=letters[nid.tip.r], d=10*nid.tip.r)
nodDt <- data.frame(c=letters[nid.int.r], e=10*nid.int.r)
## set row.names as numeric node IDs (may be changed in tests below)
row.names(allDt) <- nid.all.r
row.names(tipDt) <- nid.tip.r
row.names(nodDt) <- nid.int.r

#-----------------------------------------------------------------------
 
test.phylo4d.phylo4 <- function() {

    ## case 1: add data matching only on row position
    row.names(allDt) <- NULL
    row.names(tipDt) <- NULL
    row.names(nodDt) <- NULL

    ## these should fail because row.names don't match nodes
    checkException(phylo4d(phy.alt, tip.data=tipDt, rownamesAsLabels=TRUE))
    checkException(phylo4d(phy.alt, node.data=nodDt))

    ## brute force: no matching; with tip data
    phyd <- phylo4d(phy.alt, tip.data=tipDt, match.data=FALSE)
    checkEquals(phyd@tip.data, data.frame(tipDt,
        row.names=as.character(nid.tip)))
    checkEquals(tdata(phyd, "tip"), data.frame(tipDt,
        row.names=lab.tip))

    ## brute force: no matching; with node data
    phyd <- phylo4d(phy.alt, node.data=nodDt, match.data=FALSE)
    checkEquals(phyd@node.data, data.frame(nodDt,
        row.names=as.character(nid.int)))
    checkEquals(tdata(phyd, "internal"), data.frame(nodDt,
        row.names=lab.int))

    ## brute force: no matching; with all.data
    phyd <- phylo4d(phy.alt, all.data=allDt, match.data=FALSE)
    # TODO: these fail b/c all.data option creates numeric row.names
    # whereas tip.data and node.data options create character row.names
    #checkEquals(phyd@tip.data, data.frame(allDt,
    #    row.names=as.character(nid.all))[nid.tip,])
    #checkEquals(phyd@node.data, data.frame(allDt,
    #    row.names=as.character(nid.all))[nid.int,])
    checkEquals(tdata(phyd, "all"), data.frame(allDt,
        row.names=lab.all))

    ## brute force: no matching; with tip & node data
    ## no merging (data names don't match)
    phyd <- phylo4d(phy.alt, tip.data=tipDt["d"], node.data=nodDt["e"],
        match.data=FALSE)
    checkEquals(phyd@tip.data, data.frame(tipDt["d"], e=NA_real_,
        row.names=as.character(nid.tip)))
    checkEquals(tdata(phyd, "tip"), data.frame(tipDt["d"], e=NA_real_,
        row.names=lab.tip))
    checkEquals(phyd@node.data, data.frame(d=NA_real_, nodDt["e"],
        row.names=as.character(nid.int)))
    checkEquals(tdata(phyd, "internal"), data.frame(d=NA_real_, nodDt["e"],
        row.names=lab.int))

    ## brute force: no matching; with tip & node data
    ## merging (common data names)
    phyd <- phylo4d(phy.alt, tip.data=tipDt["c"], node.data=nodDt["c"],
        match.data=FALSE)
    checkEquals(phyd@tip.data, data.frame(c=factor(tipDt$c,
        levels=letters[nid.all]), row.names=as.character(nid.tip)))
    checkEquals(phyd@node.data, data.frame(c=factor(nodDt$c,
        levels=letters[nid.all]), row.names=as.character(nid.int)))
    checkEquals(tdata(phyd, "tip"), data.frame(c=factor(tipDt$c,
        levels=letters[nid.all]), row.names=lab.tip))
    checkEquals(tdata(phyd, "internal"), data.frame(c=factor(nodDt$c,
        levels=letters[nid.all]), row.names=lab.int))

    ## case 2: add data matching on numeric (node ID) row.names
    row.names(allDt) <- nid.all.r
    row.names(tipDt) <- nid.tip.r
    row.names(nodDt) <- nid.int.r

    ## match with node numbers, tip data
    phyd <- phylo4d(phy.alt, tip.data=tipDt)
    checkEquals(phyd@tip.data, data.frame(tipDt[order(nid.tip.r),],
        row.names=as.character(nid.tip)))
    checkEquals(tdata(phyd, "tip"), data.frame(tipDt[order(nid.tip.r),],
        row.names=lab.tip))

    ## match with node numbers, node data
    phyd <- phylo4d(phy.alt, node.data=nodDt)
    checkEquals(phyd@node.data, data.frame(nodDt[order(nid.int.r),],
        row.names=as.character(nid.int)))
    checkEquals(tdata(phyd, "internal"), data.frame(nodDt[order(nid.int.r),],
        row.names=lab.int))

    ## match with node numbers, tip & node data, no merge
    phyd <- phylo4d(phy.alt, tip.data=tipDt["d"], node.data=nodDt["e"])
    checkEquals(phyd@tip.data, data.frame(d=tipDt[order(nid.tip.r), "d"],
        e=NA_real_, row.names=as.character(nid.tip)))
    checkEquals(tdata(phyd, "tip"), data.frame(d=tipDt[order(nid.tip.r), "d"],
        e=NA_real_, row.names=lab.tip))
    checkEquals(phyd@node.data, data.frame(d=NA_real_,
        e=nodDt[order(nid.int.r), "e"], row.names=as.character(nid.int)))
    checkEquals(tdata(phyd, "internal"), data.frame(d=NA_real_,
        e=nodDt[order(nid.int.r), "e"], row.names=lab.int))

    ## match with node numbers, tip & all data
    phyd <- phylo4d(phy.alt, tip.data=tipDt, all.data=allDt)
    merged <- data.frame(merge(allDt[order(nid.all.r),],
        tipDt[order(nid.tip.r),], all=TRUE, by=0)[-1])
    checkEquals(phyd@tip.data, data.frame(merged[nid.tip,],
        row.names=as.character(nid.tip)))
    checkEquals(phyd@node.data, data.frame(merged[nid.int,],
        row.names=as.character(nid.int)))
    checkEquals(tdata(phyd, "all"), data.frame(merged, row.names=lab.all))

    ## match with node numbers, node & all data
    phyd <- phylo4d(phy.alt, node.data=nodDt, all.data=allDt)
    merged <- data.frame(merge(allDt[order(nid.all.r),],
        nodDt[order(nid.int.r),], all=TRUE, by=0)[-1])
    # TODO: need the next line because factor node data are converted
    # to character when supplied along with tip data (no merging)
    merged$c <- as.character(merged$c)
    checkEquals(phyd@tip.data, data.frame(merged[nid.tip,],
        row.names=as.character(nid.tip)))
    checkEquals(phyd@node.data, data.frame(merged[nid.int,],
        row.names=as.character(nid.int)))
    checkEquals(tdata(phyd, "all"), data.frame(merged, row.names=lab.all))

    ## match with node numbers, tip, node & all data
    phyd <- phylo4d(phy.alt, tip.data=tipDt, node.data=nodDt, all.data=allDt)
    # merge alldata with tipdata
    merged <- data.frame(merge(allDt[order(nid.all.r),],
        tipDt[order(nid.tip.r),], all=TRUE, by=0)[-1])
    # ...now merge this with nodedata
    merged <- data.frame(merge(merged, nodDt[order(nid.int.r),], all=TRUE,
        by=0, suffix=c("", ".1"))[-1])
    # TODO: need the next line because factor node data are converted
    # to character when supplied along with tip data (no merging)
    merged$c.1 <- as.character(merged$c.1)
    checkEquals(tdata(phyd, "all"), data.frame(merged, row.names=lab.all))
    # TODO: need the next line because common column names across tip
    # and node data are preserved in original form in the individual
    # data slots, but not in the tdata(x, "all") output
    names(merged)[names(merged)=="c.1"] <- "c"
    checkEquals(tdata(phyd, "tip"), data.frame(merged[nid.tip,],
        row.names=lab.tip, check.names=FALSE))
    checkEquals(tdata(phyd, "internal"), data.frame(merged[nid.int,],
        row.names=lab.int, check.names=FALSE))
    checkEquals(phyd@tip.data, data.frame(merged[nid.tip,],
        row.names=as.character(nid.tip), check.names=FALSE))
    checkEquals(phyd@node.data, data.frame(merged[nid.int,],
        row.names=as.character(nid.int), check.names=FALSE))

    ## case 3: add data matching on character (label) row.names for tips
    row.names(tipDt) <- c(lab.tip, lab.int)[nid.tip.r]
    row.names(nodDt) <- c(lab.tip, lab.int)[nid.int.r]

    ## match with names, tip data
    phyd <- phylo4d(phy.alt, tip.data=tipDt)
    checkEquals(phyd@tip.data, data.frame(tipDt[order(nid.tip.r),],
        row.names=as.character(nid.tip)))
    checkEquals(tdata(phyd, "tip"), data.frame(tipDt[order(nid.tip.r),],
        row.names=lab.tip))

    ## case 4: add data matching on mixed rowname types (for tips and
    ## for internal nodes)
    row.names(allDt)[match(nid.tip.r, nid.all.r)] <- lab.tip[nid.tip.r]
    row.names(allDt)[match(nid.int.r, nid.all.r)] <- nid.int.r

    ## match with names for tips and numbers for nodes with all data
    phyd <- phylo4d(phy.alt, all.data=allDt)
    checkEquals(tdata(phyd, "all"), data.frame(allDt[match(nid.all,
        nid.all.r),], row.names=lab.all))
    checkEquals(tdata(phyd, "tip"), data.frame(allDt[match(nid.tip,
        nid.all.r),], row.names=lab.tip))
    checkEquals(tdata(phyd, "internal"), data.frame(allDt[match(nid.int,
        nid.all.r),], row.names=lab.int))
    checkEquals(phyd@tip.data, data.frame(allDt[match(nid.tip, nid.all.r),],
        row.names=as.character(nid.tip)))
    checkEquals(phyd@node.data, data.frame(allDt[match(nid.int, nid.all.r),],
        row.names=as.character(nid.int)))

}

test.phylo4d.matrix <- function() {
}

# note: this method mostly does phylo4(phylo), then phylo4d(phylo4),
# then addData methods, which are tested more thoroughly elsewhere;
# focus here is on metadata and check.node.labels="asdata" arguments
test.phylo4d.phylo <- function() {
    # function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
    #   check.node.labels=c("keep", "drop", "asdata"), annote=list(),
    #   metadata=list(), ...)

    # show that method basically just wraps phylo4d("phylo4")
    phyd.tr <- phylo4d(tr, tip.data=tipDt, node.data=nodDt,
        all.data=allDt, match.data=TRUE, merge.data=TRUE)
    checkTrue(class(phyd.tr)=="phylo4d")
    phyd.phy <- phylo4d(phy.alt, tip.data=tipDt, node.data=nodDt,
        all.data=allDt, match.data=TRUE, merge.data=TRUE)
    # reorder for edge order consistency, then test each slot (except
    # edge labels, b/c phylo object has none)
    phyd.tr <- reorder(phyd.tr)
    phyd.phy <- reorder(phyd.phy)
    checkIdentical(edges(phyd.tr), edges(phyd.phy))
    checkIdentical(edgeLength(phyd.tr), edgeLength(phyd.phy))
    checkIdentical(nNodes(phyd.tr), nNodes(phyd.phy))
    checkIdentical(tipLabels(phyd.tr), tipLabels(phyd.phy))
    checkIdentical(nodeLabels(phyd.tr), nodeLabels(phyd.phy))
    checkIdentical(edgeOrder(phyd.tr), edgeOrder(phyd.phy))
    checkIdentical(phyd.tr@annote, phyd.phy@annote)
    # other misc checks
    checkEquals(phylo4d(phylo4(tr)), phylo4d(tr))
    checkEquals(phylo4d(phylo4(tr, check.node.labels="drop")),
        phylo4d(tr, check.node.labels="drop"))

    ##
    ## metadata
    ##

    metadata <- list(x="metadata")
    phyd <- phylo4d(tr, metadata=metadata)
    checkIdentical(metadata, phyd@metadata)

    ##
    ## check.node.labels
    ##

    # case 0: no node labels
    tr$node.label <- NULL
    phyd <- phylo4d(tr)
    checkTrue(!hasNodeLabels(phyd))

    # case 1: can't currently keep arbitrary character labels as data
    tr$node.label <- paste("n", 1:4, sep="")
    phyd <- phylo4d(tr, check.node.labels="asdata")
    checkTrue(!hasNodeLabels(phyd))
    checkTrue(all(is.na(tdata(phyd, "internal")$labelValues)))
    # the above should've produced a warning:
    opt <- options(warn=3)
    checkException(phylo4d(tr, check.node.labels="asdata"))
    options(opt)

    # case 2: convert number-like characters labels to numeric data
    tr$node.label <- as.character(1:4)
    phyd <- phylo4d(tr, check.node.labels="asdata")
    checkTrue(!hasNodeLabels(phyd))
    checkIdentical(tdata(phyd, "internal")$labelValues,
        as.numeric(tr$node.label)) 

    # case 3: convert numeric labels to numeric data
    tr$node.label <- as.numeric(1:4)
    phyd <- phylo4d(tr, check.node.labels="asdata")
    checkTrue(!hasNodeLabels(phyd))
    checkIdentical(tdata(phyd, "internal")$labelValues, tr$node.label)

    # case 4: non-unique labels can be converted to data
    tr$node.label <- rep(99, 4)
    checkException(phylo4d(tr))
    phyd <- phylo4d(tr, check.node.labels="asdata")
    checkTrue(!hasNodeLabels(phyd))
    checkIdentical(tdata(phyd, "internal")$labelValues, tr$node.label)

}

## phylo4d->phylo4d is currently unallowed
test.phylo4d.phylo4d <- function() {
    phyd <- phylo4d(phy)
    checkException(phylo4d(phyd))
}

