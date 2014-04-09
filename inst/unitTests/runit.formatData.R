#
# --- Test formatData.R ---
#

# create phylo4 object with a full complement of valid slots
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
phy.alt@label <- rev(phy@label)
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

test.formatData <- function() {
    # function(phy, dt, type=c("tip", "internal", "all"),
    #   match.data=TRUE, rownamesAsLabels=FALSE,
    #   label.type=c("rownames", "column"), label.column=1,
    #   missing.data=c("fail", "warn", "OK"),
    #   extra.data=c("warn", "OK", "fail"), keep.all=TRUE

    ## vector data coerced to data.frame (colname dt)
    checkIdentical(phylobase:::formatData(phy.alt, 1:5),
        phylobase:::formatData(phy.alt, data.frame(dt=1:5)))
    ## list of vector data coerced to data.frame (colnames as given)
    checkIdentical(phylobase:::formatData(phy.alt, list(a=1:5, b=6:10)),
        phylobase:::formatData(phy.alt, data.frame(a=1:5, b=6:10)))
    ## factor data coerced to data.frame (colname dt)
    checkIdentical(phylobase:::formatData(phy.alt, factor(letters[1:5])),
        phylobase:::formatData(phy.alt, data.frame(dt=letters[1:5])))
    ## matrix data coerced to data.frame (colnames V1, V2)
    checkIdentical(phylobase:::formatData(phy.alt, matrix(1:10, ncol=2)),
        phylobase:::formatData(phy.alt, data.frame(V1=1:5, V2=6:10)))
    ## matrix data coerced to data.frame (colname as given)
    checkIdentical(phylobase:::formatData(phy.alt, matrix(1:10, ncol=2,
        dimnames=list(NULL, c("a", "b")))),
        phylobase:::formatData(phy.alt, data.frame(a=1:5, b=6:10)))
    ## error if dt is, say, a phylo4 object
    checkException(phylobase:::formatData(phy.alt, phy.alt))

    ## error if column number is out of range
    checkException(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=FALSE,
        label.type="column", label.column=3),
        data.frame(a=c(1:5, rep(NA, 4)), lab=c(rev(nid.tip), rep(NA,
        4)), row.names=nid.all))
    ## error if column name is wrong
    checkException(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=FALSE,
        label.type="column", label.column="foo"),
        data.frame(a=c(1:5, rep(NA, 4)), lab=c(rev(nid.tip), rep(NA,
        4)), row.names=nid.all))


    #
    # matching options
    #

    ## don't match (purely positional)
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip", match.data=FALSE),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    ## match on rownames (node numbers)
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip", match.data=TRUE),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip"), data.frame(a=c(5:1,
        rep(NA, 4)), row.names=nid.all))
    ## match on rownames (labels)
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=rev(lab.tip)), type="tip", match.data=TRUE),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## match on rownames (mixed node numbers and labels)
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## but fails if rownamesAsLabels is TRUE
    checkException(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE, rownamesAsLabels=TRUE))

    #
    #   label.type="column" and label.column=2
    #

    ## should ignore label (purely positional) and retain a label col
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=FALSE,
        label.type="column", label.column=2),
        data.frame(a=c(1:5, rep(NA, 4)), lab=c(rev(nid.tip), rep(NA,
        4)), row.names=nid.all))
    ## match on label column (node numbers)
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=TRUE,
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip",
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## match on label column (labels)
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(lab.tip)), type="tip", match.data=TRUE,
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(lab.tip)), type="tip", match.data=TRUE,
        label.type="column", label.column="lab"),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## match on label column (mixed node numbers and labels)
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])), type="tip",
        match.data=TRUE, label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## but fails if rownamesAsLabels is TRUE
    checkException(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE, rownamesAsLabels=TRUE,
        label.type="column", label.column=2))

    ## try to match internal nodes when type='tips'
    checkException(phylobase:::formatData(phy.alt, data.frame(a=1:5, row.names=4:8),
        type="tip"))
    ## and vice versa
    checkException(phylobase:::formatData(phy.alt, data.frame(a=6:9, row.names=1:4),
        type="internal"))

    #
    # missing.data
    #

    ## force error conditions
    checkException(phylobase:::formatData(phy.alt, data.frame(a=1:3), type="tip"))
    checkException(phylobase:::formatData(phy.alt, data.frame(a=1:3), type="tip",
        missing.data="fail"))
    options(warn=3)
    checkException(phylobase:::formatData(phy.alt, data.frame(a=1:3), type="tip",
        missing.data="warn"))
    options(warn=0)
    ## missing data with matching
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.tip)[-1],
        row.names=rev(nid.tip)[-1]), type="tip", missing.data="OK"),
        data.frame(a=c(nid.tip[-5], rep(NA, 5))))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.int)[-1],
        row.names=rev(nid.int)[-1]), type="internal", missing.data="OK"),
        data.frame(a=c(rep(NA, 5), nid.int[-4], NA)))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.all)[-1],
        row.names=rev(nid.all)[-1]), type="all", missing.data="OK"),
        data.frame(a=c(nid.all[-9], NA)))
    ## missing data without matching
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.tip)[-1]),
        type="tip", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rev(nid.tip)[-1], rep(NA, 5))))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.int)[-1]),
        type="internal", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rep(NA, 5), rev(nid.int)[-1], NA)))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.all)[-1]),
        type="all", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rev(nid.all)[-1], NA)))

    #
    # extra.data
    #

    ## force error conditions
    checkException(phylobase:::formatData(phy.alt, data.frame(a=1:3), type="tip",
        missing.data="fail"))
    options(warn=3)
    checkException(phylobase:::formatData(phy.alt, data.frame(a=0:5, row.names=0:5),
        type="tip", missing="warn"))
    checkException(phylobase:::formatData(phy.alt, data.frame(a=0:5, row.names=0:5),
        type="tip"))
    options(warn=0)
    ## extra data with matching
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=c(0L, rev(nid.tip)),
        row.names=c(0, rev(nid.tip))), type="tip", extra.data="OK"),
        data.frame(a=c(nid.tip, rep(NA, 4))))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=c(0L, rev(nid.int)),
        row.names=c(0, rev(nid.int))), type="internal", extra.data="OK"),
        data.frame(a=c(rep(NA, 5), nid.int)))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=c(0L, rev(nid.all)),
        row.names=c(0, rev(nid.all))), type="all", extra.data="OK"),
        data.frame(a=nid.all))
    ## extra data without matching
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:15),
        type="tip", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(1:5, rep(NA, 4))))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:15),
        type="internal", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(rep(NA, 5), 1:4)))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:15),
        type="all", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(1:9)))

    ## allow both extra.data and missing.data
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=0:3, row.names=0:3),
        type="tip", extra.data="OK", missing.data="OK"),
        data.frame(a=c(1:3, rep(NA, 6))))

    #
    # keep.all
    #

    ## keep all rows
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=nid.tip), type="tip", keep.all=TRUE),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=nid.tip), type="tip"),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=6:9,
        row.names=nid.int), type="internal", keep.all=TRUE),
        data.frame(a=c(rep(NA, 5), 6:9), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=6:9,
        row.names=nid.int), type="internal"),
        data.frame(a=c(rep(NA, 5), 6:9), row.names=nid.all))
    ## only keep 'type' rows
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=nid.tip), type="tip", keep.all=FALSE),
        data.frame(a=c(1:5), row.names=nid.tip))
    checkIdentical(phylobase:::formatData(phy.alt, data.frame(a=6:9,
        row.names=nid.int), type="internal", keep.all=FALSE),
        data.frame(a=c(6:9), row.names=nid.int))

}

test.formatDataWithDup <- function() {

    ## Saving default options
    op <- phylobase.options()

    ## Changing default options
    phylobase.options(allow.duplicated.labels="ok") 

    ## Creating phylo4 object with duplicated labels
    phy.dup <- phy.alt
    tipLabels(phy.dup)[2] <- tipLabels(phy.dup)[1]

    ## vector data coerced to data.frame (colname dt)
    checkIdentical(phylobase:::formatData(phy.dup, 1:5),
        phylobase:::formatData(phy.dup, data.frame(dt=1:5)))
    ## list of vector data coerced to data.frame (colnames as given)
    checkIdentical(phylobase:::formatData(phy.dup, list(a=1:5, b=6:10)),
        phylobase:::formatData(phy.dup, data.frame(a=1:5, b=6:10)))
    ## factor data coerced to data.frame (colname dt)
    checkIdentical(phylobase:::formatData(phy.dup, factor(letters[1:5])),
        phylobase:::formatData(phy.dup, data.frame(dt=letters[1:5])))
    ## matrix data coerced to data.frame (colnames V1, V2)
    checkIdentical(phylobase:::formatData(phy.dup, matrix(1:10, ncol=2)),
        phylobase:::formatData(phy.dup, data.frame(V1=1:5, V2=6:10)))
    ## matrix data coerced to data.frame (colname as given)
    checkIdentical(phylobase:::formatData(phy.dup, matrix(1:10, ncol=2,
        dimnames=list(NULL, c("a", "b")))),
        phylobase:::formatData(phy.dup, data.frame(a=1:5, b=6:10)))
    ## error if dt is, say, a phylo4 object
    checkException(phylobase:::formatData(phy.dup, phy.dup))

    #
    # matching options
    #

    ## don't match (purely positional)
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip", match.data=FALSE),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    ## match on rownames (node numbers)
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip", match.data=TRUE),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip"), data.frame(a=c(5:1,
        rep(NA, 4)), row.names=nid.all))
    ## match on rownames (labels)
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=c(1,3,4,5),
        row.names=rev(lab.tip[-2])), type="tip", match.data=TRUE),
        data.frame(a=c(5,5,4,3,1, rep(NA, 4)), row.names=nid.all))
    ## match on rownames (mixed node numbers and labels)
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=c(1,2,3,4,5),
        row.names=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE),
        data.frame(a=c(5,4,3,2,1, rep(NA, 4)), row.names=nid.all))
    ## but fails if rownamesAsLabels is TRUE
    checkException(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE, rownamesAsLabels=TRUE))

    #
    #   label.type="column" and label.column=2
    #

    ## should ignore label (purely positional) and retain a label col
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=FALSE,
        label.type="column", label.column=2),
        data.frame(a=c(1:5, rep(NA, 4)), lab=c(rev(nid.tip), rep(NA,
        4)), row.names=nid.all))
    ## match on label column (node numbers)
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=TRUE,
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip",
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## match on label column (labels)
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:4,
        lab=rev(lab.tip[-2])), type="tip", match.data=TRUE,
        label.type="column", label.column=2),
        data.frame(a=as.integer(c(4, 4:1, rep(NA, 4))), row.names=nid.all))
    ## match on label column (mixed node numbers and labels)
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])), type="tip",
        match.data=TRUE, label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## but fails if rownamesAsLabels is TRUE
    checkException(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE, rownamesAsLabels=TRUE,
        label.type="column", label.column=2))

    ## try to match internal nodes when type='tips'
    checkException(phylobase:::formatData(phy.dup, data.frame(a=1:5, row.names=4:8),
        type="tip"))
    ## and vice versa
    checkException(phylobase:::formatData(phy.dup, data.frame(a=6:9, row.names=1:4),
        type="internal"))

    #
    # missing.data
    #

    ## force error conditions
    checkException(phylobase:::formatData(phy.dup, data.frame(a=1:3), type="tip"))
    checkException(phylobase:::formatData(phy.dup, data.frame(a=1:3), type="tip",
        missing.data="fail"))
    options(warn=3)
    checkException(phylobase:::formatData(phy.dup, data.frame(a=1:3), type="tip",
        missing.data="warn"))
    options(warn=0)
    ## missing data with matching
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.tip)[-1],
        row.names=rev(nid.tip)[-1]), type="tip", missing.data="OK"),
        data.frame(a=c(nid.tip[-5], rep(NA, 5))))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.int)[-1],
        row.names=rev(nid.int)[-1]), type="internal", missing.data="OK"),
        data.frame(a=c(rep(NA, 5), nid.int[-4], NA)))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.all)[-1],
        row.names=rev(nid.all)[-1]), type="all", missing.data="OK"),
        data.frame(a=c(nid.all[-9], NA)))
    ## missing data without matching
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.tip)[-1]),
        type="tip", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rev(nid.tip)[-1], rep(NA, 5))))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.int)[-1]),
        type="internal", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rep(NA, 5), rev(nid.int)[-1], NA)))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.all)[-1]),
        type="all", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rev(nid.all)[-1], NA)))

    #
    # extra.data
    #

   ## force error conditions
    checkException(phylobase:::formatData(phy.dup, data.frame(a=1:3), type="tip",
        missing.data="fail"))
    options(warn=3)
    checkException(phylobase:::formatData(phy.dup, data.frame(a=0:5, row.names=0:5),
        type="tip", missing="warn"))
    checkException(phylobase:::formatData(phy.dup, data.frame(a=0:5, row.names=0:5),
        type="tip"))
    options(warn=0)
    ## extra data with matching
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=c(0L, rev(nid.tip)),
        row.names=c(0, rev(nid.tip))), type="tip", extra.data="OK"),
        data.frame(a=c(nid.tip, rep(NA, 4))))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=c(0L, rev(nid.int)),
        row.names=c(0, rev(nid.int))), type="internal", extra.data="OK"),
        data.frame(a=c(rep(NA, 5), nid.int)))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=c(0L, rev(nid.all)),
        row.names=c(0, rev(nid.all))), type="all", extra.data="OK"),
        data.frame(a=nid.all))
    ## extra data without matching
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:15),
        type="tip", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(1:5, rep(NA, 4))))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:15),
        type="internal", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(rep(NA, 5), 1:4)))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:15),
        type="all", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(1:9)))

    ## allow both extra.data and missing.data
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=0:3, row.names=0:3),
        type="tip", extra.data="OK", missing.data="OK"),
        data.frame(a=c(1:3, rep(NA, 6))))

    #
    # keep.all
    #

    ## keep all rows
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=nid.tip), type="tip", keep.all=TRUE),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=nid.tip), type="tip"),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=6:9,
        row.names=nid.int), type="internal", keep.all=TRUE),
        data.frame(a=c(rep(NA, 5), 6:9), row.names=nid.all))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=6:9,
        row.names=nid.int), type="internal"),
        data.frame(a=c(rep(NA, 5), 6:9), row.names=nid.all))
    ## only keep 'type' rows
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=nid.tip), type="tip", keep.all=FALSE),
        data.frame(a=c(1:5), row.names=nid.tip))
    checkIdentical(phylobase:::formatData(phy.dup, data.frame(a=6:9,
        row.names=nid.int), type="internal", keep.all=FALSE),
        data.frame(a=c(6:9), row.names=nid.int))

    ## restoring default options
    phylobase.options(op)
}
