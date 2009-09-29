data(geospiza)

p4 <- extractTree(geospiza)
p4d <- geospiza

test.labelsTipsPhylo4 <- function() {
    tLbl <- paste("t", 1:nTips(p4), sep="")
    nmTLbl <- tLbl
    names(nmTLbl) <- nodeId(p4, "tip")
    tLbl <- sample(tLbl)
    nmTLbl <- sample(nmTLbl)

    ## case all options by default and unnamed vector
    p4c <- p4
    tipLabels(p4c) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4c)) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(tipLabels(p4c)), tLbl)

    ## case all options by default and named vector
    p4c <- p4
    tipLabels(p4c) <- nmTLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4c)) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(tipLabels(p4c)), unname(nmTLbl))

    ## case type defined
    p4c <- p4
    tipLabels(p4c) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4c)) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(tipLabels(p4c)), tLbl)

    ## case type defined and use.names=TRUE but no names
    p4c <- p4
    tipLabels(p4c, use.names=TRUE) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4c)) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(tipLabels(p4c)), tLbl)

    ## case type defined and use.names=TRUE with names
    p4c <- p4
    tipLabels(p4c, use.names=TRUE) <- nmTLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4c)) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(tipLabels(p4c), nmTLbl[order(as.numeric(names(nmTLbl)))])
}

test.labelsNodePhylo4 <- function() {

    ndLbl <- paste("n", 1:nNodes(p4), sep="")
    nmNdLbl <- ndLbl
    names(nmNdLbl) <- nodeId(p4, "internal")

    ndLbl <- sample(ndLbl)
    nmNdLbl <- sample(nmNdLbl)

    ## case type defined
    p4c <- p4
    nodeLabels(p4c) <- ndLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(nodeLabels(p4c)) %in% nodeId(p4c, "all")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(nodeLabels(p4c)), ndLbl)

    ## case type defined and use.names=TRUE but no names
    p4c <- p4
    nodeLabels(p4c, use.names=TRUE) <- ndLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(nodeLabels(p4c)) %in% nodeId(p4c, "all")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(nodeLabels(p4c)), ndLbl)

    ## case type defined and use.names=TRUE with names
    p4c <- p4
    nodeLabels(p4c, use.names=TRUE) <- nmNdLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(nodeLabels(p4c)) %in% nodeId(p4c, "all")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(nodeLabels(p4c), nmNdLbl[order(as.numeric(names(nmNdLbl)))])

}

test.labelsTipsPhylo4d <- function() {
    tLbl <- paste("t", 1:nTips(p4d), sep="")
    nmTLbl <- tLbl
    names(nmTLbl) <- nodeId(p4d, "tip")
    tLbl <- sample(tLbl)
    nmTLbl <- sample(nmTLbl)

    ## case all options by default and unnamed vector
    p4dc <- p4d
    tipLabels(p4dc) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4dc)) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(tipLabels(p4dc)), tLbl)

    ## case all options by default and named vector
    p4dc <- p4d
    tipLabels(p4dc) <- nmTLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4dc)) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(tipLabels(p4dc)), unname(nmTLbl))

    ## case type defined
    p4dc <- p4d
    tipLabels(p4dc) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4dc)) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(tipLabels(p4dc)), tLbl)

    ## case type defined and use.names=TRUE but no names
    p4dc <- p4d
    tipLabels(p4dc, use.names=TRUE) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4dc)) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(tipLabels(p4dc)), tLbl)

    ## case type defined and use.names=TRUE with names
    p4dc <- p4d
    tipLabels(p4dc, use.names=TRUE) <- nmTLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(tipLabels(p4dc)) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(tipLabels(p4dc), nmTLbl[order(as.numeric(names(nmTLbl)))])
}

test.labelsNodePhylo4d <- function() {

    ndLbl <- paste("n", 1:nNodes(p4d), sep="")
    nmNdLbl <- ndLbl
    names(nmNdLbl) <- nodeId(p4d, "internal")

    ndLbl <- sample(ndLbl)
    nmNdLbl <- sample(nmNdLbl)

    ## case type defined
    p4dc <- p4d
    nodeLabels(p4dc) <- ndLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(nodeLabels(p4dc)) %in% nodeId(p4dc, "internal")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(nodeLabels(p4dc)), ndLbl)

    ## case type defined and use.names=TRUE but no names
    p4dc <- p4d
    nodeLabels(p4dc, use.names=TRUE) <- ndLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(nodeLabels(p4dc)) %in% nodeId(p4dc, "internal")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(nodeLabels(p4dc)), ndLbl)

    ## case type defined and use.names=TRUE with names
    p4dc <- p4d
    nodeLabels(p4dc, use.names=TRUE) <- nmNdLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(nodeLabels(p4dc)) %in% nodeId(p4dc, "internal")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(nodeLabels(p4dc), nmNdLbl[order(as.numeric(names(nmNdLbl)))])

}

test.labelsAllPhylo4 <- function() {

    allLbl <- paste("n", 1:(nTips(p4)+nNodes(p4)), sep="")
    nmAllLbl <- allLbl
    names(nmAllLbl) <- nodeId(p4, "all")

    allLbl <- sample(allLbl)
    nmAllLbl <- sample(nmAllLbl)

    p4c <- p4
    labels(p4c, "all") <- allLbl
    checkTrue(all(names(labels(p4c)) %in% nodeId(p4c, "all")))
    checkEquals(unname(labels(p4c)), allLbl)

    p4c <- p4
    labels(p4c, "all") <- nmAllLbl
    checkTrue(all(names(labels(p4c)) %in% nodeId(p4c, "all")))
    checkEquals(unname(labels(p4c)), unname(nmAllLbl))

    p4c <- p4
    tmpNm <- nmAllLbl[order(as.numeric(names(nmAllLbl)))]
    labels(p4c, "all", use.names=TRUE) <- nmAllLbl
    checkTrue(all(names(labels(p4c)) %in% nodeId(p4c, "all")))
    checkEquals(labels(p4c), tmpNm[names(tmpNm) %in% nodeId(p4c, "all")])
}

test.labelsAllPhylo4d <- function() {

    allLbl <- paste("n", 1:(nTips(p4d)+nNodes(p4d)), sep="")
    nmAllLbl <- allLbl
    names(nmAllLbl) <- nodeId(p4d, "all")

    allLbl <- sample(allLbl)
    nmAllLbl <- sample(nmAllLbl)

    p4dc <- p4d
    labels(p4dc, "all") <- allLbl
    checkTrue(all(names(labels(p4dc)) %in% nodeId(p4dc, "all")))
    checkEquals(unname(labels(p4dc)), allLbl)

    p4dc <- p4d
    labels(p4dc, "all") <- nmAllLbl
    checkTrue(all(names(labels(p4dc)) %in% nodeId(p4dc, "all")))
    checkEquals(unname(labels(p4dc)), unname(nmAllLbl))


    p4dc <- p4d
    tmpNm <- nmAllLbl[order(as.numeric(names(nmAllLbl)))]
    labels(p4dc, "all", use.names=TRUE) <- nmAllLbl
    checkTrue(all(names(labels(p4dc)) %in% nodeId(p4dc, "all")))
    checkEquals(labels(p4dc), tmpNm[names(tmpNm) %in% nodeId(p4dc, "all")])
}
