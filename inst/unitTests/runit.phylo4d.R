test.Phylo4d.bruteforce <- function() {
    data(geospiza_raw)
    ## using the raw data
    tr <- geospiza_raw$tree
    tr <- as(tr, "phylo4")

    ## creating some data with rownames matching tip numbers
    tDt <- data.frame(testData1 = 1:nTips(tr),
                      testData2 = (1:nTips(tr))*2)
    tDt <- tDt[sample(1:nrow(tDt)),]

    nDt <- data.frame(testData1 = 1:nNodes(tr))
    nDt <- nDt[sample(1:nrow(nDt)),]

    aDt <- data.frame(testData1 = 1:(nNodes(tr)+nTips(tr)))
    aDt <- aDt[sample(1:nrow(aDt)) ,, drop=FALSE]

    ## brute force: no matching; with tip data
    xx <- phylo4d(tr, tip.data=tDt, match.data=FALSE)
    checkEquals(xx@tip.data[,1], tDt[,1])
    checkEquals(tdata(xx)[,1], tDt[,1])

    ## brute force: no matching; with node data
    yy <- phylo4d(tr, node.data=nDt, match.data=FALSE)
    checkEquals(yy@node.data[,1], nDt)
    checkEquals(tdata(yy, "internal")[,1], nDt)

    ## brute force: no matching; with all.data
    zz <- phylo4d(tr, all.data=aDt, match.data=FALSE)
    checkEquals(zz@tip.data[,1], aDt[nodeId(tr, "tip"),1])
    checkEquals(zz@node.data[,1], aDt[nodeId(tr, "internal"),1])
    checkEquals(tdata(zz, "all")[,1], aDt[,1])

    ## brute force: no matching; with tip & node data
    ## no merging (data names don't match)
    xx <- phylo4d(tr, tip.data=tDt, node.data=nDt, match.data=FALSE)
    checkEquals(xx@tip.data[,1], tDt[,1])
    checkEquals(xx@node.data[,3], nDt)
    checkEquals(tdata(xx, "tip")[,1], tDt[,1])
    checkEquals(tdata(xx, "internal")[,3], nDt)

    ## brute force: no matching; with tip & node data
    ## merging
    nDt <- data.frame(nDt)
    names(nDt) <- names(tDt)[1]
    xx <- phylo4d(tr, tip.data=tDt[,1,drop=F], node.data=nDt, match.data=FALSE)
    checkEquals(xx@tip.data[,1], tDt[,1])
    checkEquals(xx@node.data[,1], nDt[,1])
    checkEquals(tdata(xx, "tip")[,1], tDt[,1])
    checkEquals(tdata(xx, "internal")[,1], nDt[,1])

}

test.Phylo4d.withNb <- function() {
    data(geospiza_raw)
    ## using the raw data
    tr <- geospiza_raw$tree
    tr <- as(tr, "phylo4")

    ## creating some data with rownames matching tip numbers
    tDt <- data.frame(testData1 = 1:nTips(tr),
                      testData2 = (1:nTips(tr))*2)
    tDt <- tDt[sample(1:nrow(tDt)),]

    nDt <- data.frame(testData1 = 1:nNodes(tr))
    rownames(nDt) <- nodeId(tr, "internal")
    nDt <- nDt[sample(1:nrow(nDt)) ,, drop=FALSE]

    aDt <- data.frame(testData1 = 1:(nNodes(tr)+nTips(tr)))
    aDt <- aDt[sample(1:nrow(aDt)) ,, drop=FALSE]

    ## match with node numbers, tip data
    xx <- phylo4d(tr, tip.data=tDt)
    checkEquals(xx@tip.data[,1], 1:nTips(tr))
    checkEquals(tdata(xx, "tip")[,1], 1:nTips(tr))

    ## match with node numbers, node data
    xx <- phylo4d(tr, node.data=nDt)
    checkEquals(xx@node.data[,1], 1:nNodes(tr))
    checkEquals(tdata(xx, "internal")[,1], 1:nNodes(tr))

    ## match with node numbers, tip & node data
    xx <- phylo4d(tr, tip.data=tDt, node.data=nDt)
    checkEquals(xx@tip.data[,1], 1:nTips(tr))
    checkEquals(xx@node.data[,3], 1:nNodes(tr))
    checkEquals(tdata(xx, "tip")[,1], 1:nTips(tr))
    checkEquals(tdata(xx, "internal")[,3], 1:nNodes(tr))

    ## match with node numbers, tip & all data
    xx <- phylo4d(tr, tip.data=tDt, all.data=aDt)
    checkEquals(xx@tip.data[,1], 1:nTips(tr))
    checkEquals(xx@node.data[,1], (nTips(tr)+1):(nTips(tr)+nNodes(tr)))
    checkEquals(tdata(xx, "all")[,1], 1:(nTips(tr)+nNodes(tr)))

    ## match with node numbers, node & all data
    xx <- phylo4d(tr, node.data=nDt, all.data=aDt)
    checkEquals(xx@tip.data[,1], 1:nTips(tr))
    checkEquals(xx@node.data[,2], 1:nNodes(tr))
    checkEquals(tdata(xx, "all")[,1], 1:(nTips(tr)+nNodes(tr)))

    ## match with node numbers, tip, node & all data
    xx <- phylo4d(tr, tip.data=tDt, node.data=nDt, all.data=aDt)
    checkEquals(xx@tip.data[,1], 1:nTips(tr))
    checkEquals(xx@tip.data[,2], 1:nTips(tr))
    checkEquals(xx@node.data[,1], (nTips(tr)+1):(nTips(tr)+nNodes(tr)))
    checkEquals(xx@node.data[,4], 1:nNodes(tr))
    checkEquals(tdata(xx, "all")[,1], 1:(nTips(tr)+nNodes(tr)))
}

test.Phylo4d.withNames <- function() {
    data(geospiza_raw)
    ## using the raw data
    tr <- geospiza_raw$tree
    tr <- as(tr, "phylo4")

    ## creating some data with rownames matching tip numbers
    tDt <- data.frame(testData1 = 1:nTips(tr),
                      testData2 = (1:nTips(tr))*2)
    rownames(tDt) <- tipLabels(tr)
    tDt <- tDt[sample(1:nrow(tDt)),]

    aDt <- data.frame(testData1 = 1:(nNodes(tr)+nTips(tr)))
    rownames(aDt)[1:nTips(tr)] <- tipLabels(tr)
    aDt <- aDt[sample(1:nrow(aDt)) ,, drop=FALSE]

    ## match with names, tip data
    xx <- phylo4d(tr, tip.data=tDt)
    checkEquals(xx@tip.data[,1], 1:nTips(tr))
    checkEquals(tdata(xx, "tip")[,1], 1:nTips(tr))

    ## match with names for tips and numbers for nodes with all data
    xx <- phylo4d(tr, all.data=aDt)
    checkEquals(xx@tip.data[,1], 1:nTips(tr))
    checkEquals(xx@node.data[,1], (nTips(tr)+1):(nTips(tr)+nNodes(tr)))
    checkEquals(tdata(xx, "all")[,1], 1:(nTips(tr)+nNodes(tr)))
}
