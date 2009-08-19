require(phylobase)
data(geospiza)

gtree <- extractTree(geospiza)
stopifnot(identical(gtree,prune(gtree,character(0))))

stopifnot(identical(tdata(subset(geospiza)),
                    tdata(subset(geospiza, tipLabels(geospiza)))))



