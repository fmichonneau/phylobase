#
# --- Test readNexus.R ---
#

### Get all the test files
if (Sys.getenv("RCMDCHECK") == FALSE) {
    pth <- file.path(getwd(), "..", "inst", "nexusfiles")
} else {
    pth <- system.file(package="phylobase", "nexusfiles")
}
## co1.nex -- typical output from MrBayes. Contains 2 identical trees, the first
## one having posterior probabilities as node labels
co1File <- file.path(pth, "co1.nex")

## MultiLineTrees.nex -- 2 identical trees stored on several lines
multiLinesFile <- file.path(pth, "MultiLineTrees.nex")

## treepluscharV01.nex -- Mesquite file
treepluschar <- file.path(pth, "treepluscharV01.nex")

## Contain correct (as of 2010-03-08) phylo4 representation of one of the tree
## stored in the nexus file
mlFile <- file.path(pth, "multiLines.Rdata")

stopifnot(file.exists(co1File))
stopifnot(file.exists(treepluschar))
stopifnot(file.exists(multiLinesFile))
stopifnot(file.exists(mlFile))

test.readNexus <- function() {
    ## function (file, simplify=TRUE, type=c("all", "tree", "data"),
    ##   char.all=FALSE, polymorphic.convert=TRUE, levels.uniform=TRUE,
    ##   check.node.labels=c("keep", "drop", "asdata"))

    ## ########### CO1 -- MrBayes file
    ## Tree properties
    ## Labels
    labCo1 <- c("Cow", "Seal", "Carp", "Loach", "Frog", "Chicken", "Human",
                "Mouse", "Rat", "Whale", NA, NA, NA, NA, NA, NA, NA, NA)
    names(labCo1) <- 1:18
    ## Edge lengths
    eLco1 <- c(0.143336, 0.225087, 0.047441, 0.055934, 0.124549, 0.204809, 0.073060, 0.194575,
               0.171296, 0.222039, 0.237101, 0.546258, 0.533183, 0.154442, 0.134574, 0.113163,
               0.145592)
    names(eLco1) <- c("11-1", "11-2", "11-12", "12-13", "13-14", "14-15", "15-16", "16-17", "17-3",
                      "17-4", "16-5", "15-6", "14-7", "13-18", "18-8", "18-9", "12-10")
    ## Node types
    nTco1 <-  c("tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip",
                "tip", "internal", "internal", "internal", "internal", "internal",
                "internal", "internal", "internal")
    names(nTco1) <- 1:18
    ## Label values
    lVco1 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.93, 0.88, 0.99, 1.00,
               0.76, 1.00, 1.00)
    ## Read trees
    co1 <- readNexus(file=co1File, check.node.labels="asdata")
    ## Tree 1
    co1Tree1 <- co1[[1]]
    checkIdentical(labels(co1Tree1), labCo1)     # check labels
    checkIdentical(edgeLength(co1Tree1), eLco1)  # check edge lengths
    checkIdentical(nodeType(co1Tree1), nTco1)    # check node types
    checkIdentical(as(co1Tree1, "data.frame")$labelValues, lVco1) # check label values
    ## Tree 2
    co1Tree2 <- co1[[2]]
    checkIdentical(labels(co1Tree2), labCo1)     # check labels
    checkIdentical(edgeLength(co1Tree2), eLco1)  # check edge lengths
    checkIdentical(nodeType(co1Tree2), nTco1)    # check node types

    ## Check option simplify
    co1 <- readNexus(file=co1File, check.node.labels="asdata", simplify=TRUE)
    checkIdentical(length(co1), as.integer(1))   # make sure there is only one tree
    checkIdentical(labels(co1), labCo1)          # check labels
    checkIdentical(edgeLength(co1), eLco1)       # check edge lengths
    checkIdentical(nodeType(co1), nTco1)         # check node type
    checkIdentical(as(co1, "data.frame")$labelValues, lVco1)  # check label values

    ## Check option check.node.labels
    checkException(readNexus(file=co1File, check.node.labels="keep")) # fail because labels aren't unique
    co1 <- readNexus(file=co1File, check.node.labels="drop", simplify=TRUE)
    checkIdentical(labels(co1), labCo1)          # check labels
    checkIdentical(edgeLength(co1), eLco1)       # check edge lengths
    checkIdentical(nodeType(co1), nTco1)         # check node type
    checkIdentical(as(co1, "data.frame")$labelValues, NULL)  # check label values don't exist

    ## ########### Mutli Lines
    multiLines <- readNexus(file=multiLinesFile)
    ## load correct representation and make sure that the trees read
    ## match it
    load(mlFile)
    checkIdentical(multiLines[[1]], ml1)
    checkIdentical(multiLines[[2]], ml1)
    rm(ml1)

    ## ########### Tree + data -- file from Mesquite
    ## tree properties
    labTr <-  c("Myrmecocystussemirufus", "Myrmecocystusplacodops",
                "Myrmecocystusmendax", "Myrmecocystuskathjuli",
                "Myrmecocystuswheeleri", "Myrmecocystusmimicus",
                "Myrmecocystusdepilis", "Myrmecocystusromainei",
                "Myrmecocystusnequazcatl", "Myrmecocystusyuma",
                "Myrmecocystuskennedyi", "Myrmecocystuscreightoni",
                "Myrmecocystussnellingi", "Myrmecocystustenuinodis",
                "Myrmecocystustestaceus", "Myrmecocystusmexicanus",
                "Myrmecocystuscfnavajo", "Myrmecocystusnavajo",
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    names(labTr) <- 1:35
    eTr <- c(1.699299, 0.894820, 0.836689, 4.524387, 0.506099, 0.198842, 0.689044,
             2.926053, 1.724765, 1.724765, 4.650818, 4.255993, 1.083870, 1.083870,
             0.802512, 2.027251, 2.708942, 2.708942, NA, 0.284767, 2.257581,
             2.193845, 2.193845, 4.451425, 6.044804, 10.569191, 8.635503, 2.770378,
             2.770378, 12.300701, 8.275077, 5.724923, 2.855375, 2.869547, 2.869547)
    names(eTr) <- c("19-20","20-21","21-22","22-23","23-24","24-25","25-26","26-27",
                    "27-1", "27-2","26-3","25-28","28-4","28-5","24-29","29-30",
                    "30-6","30-7","0-19","29-31","31-32","32-8","32-9","31-10",
                    "23-11","22-12","21-33","33-13","33-14","20-15","19-34","34-16",
                    "34-35","35-17","35-18")
    nTtr <- c("tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip",
              "tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip",
              "root", "internal", "internal", "internal", "internal", "internal",
              "internal", "internal", "internal", "internal", "internal",
              "internal", "internal", "internal", "internal", "internal",
              "internal")
    names(nTtr) <- 1:35
    ## data to test against
    dtTest1 <- data.frame(time = factor(c(2,1,0,0,0,0,2,0,2,0,0,0,0,1,1,1,0,1)),
                          subgenus = factor(c(2,1,0,0,0,0,2,0,2,0,0,0,0,1,1,2,0,1)))
    row.names(dtTest1) <- c("Myrmecocystuscfnavajo","Myrmecocystuscreightoni",
                            "Myrmecocystusdepilis","Myrmecocystuskathjuli",
                            "Myrmecocystuskennedyi","Myrmecocystusmendax",
                            "Myrmecocystusmexicanus","Myrmecocystusmimicus",
                            "Myrmecocystusnavajo","Myrmecocystusnequazcatl",
                            "Myrmecocystusplacodops","Myrmecocystusromainei",
                            "Myrmecocystussemirufus","Myrmecocystussnellingi",
                            "Myrmecocystustenuinodis","Myrmecocystustestaceus",
                            "Myrmecocystuswheeleri","Myrmecocystusyuma")
    dtTest2 <- dtTest1
    dtTest2$time <- ifelse(dtTest2$time == 0, "diurnal",
                           ifelse(dtTest2$time == 1, "crepuscular", "nocturnal"))
    dtTest2$time <- factor(dtTest2$time,
                           levels=c("diurnal", "crepuscular", "nocturnal"))
    dtTest2$subgenus <- ifelse(dtTest2$subgenus == 0, "Endiodioctes",
                               ifelse(dtTest2$subgenus == 1, "Eremnocystus", "Myrmecocystus"))
    dtTest2$subgenus <- factor(dtTest2$subgenus)
    p4 <- "phylo4"
    p4d <- "phylo4d"
    attributes(p4) <- attributes(p4d) <- list(package="phylobase")
    ## Tree only
    tr <- readNexus(file=treepluschar, type="tree")
    checkIdentical(labels(tr), labTr)   # check labels
    checkIdentical(edgeLength(tr), eTr) # check edge lengths
    checkIdentical(nodeType(tr), nTtr)  # check node types
    checkIdentical(class(tr), p4)       # check class
    ## Data only
    dt1 <- readNexus(file=treepluschar, type="data")
    checkIdentical(dt1, dtTest1)
    dt2 <- readNexus(file=treepluschar, type="data", levels.uniform=FALSE)
    checkIdentical(dt2, dtTest2)
    ## Tree + Data
    trDt1 <- readNexus(file=treepluschar, type="all")
    str(trDt1)
    print(labels(trDt1))
    print(labTr)
    checkIdentical(labels(trDt1), labTr)   # check labels
    checkIdentical(edgeLength(trDt1), eTr) # check edge lengths
    checkIdentical(nodeType(trDt1), nTtr)  # check node types
    checkIdentical(class(trDt1), p4d)      # check class
    checkIdentical(tdata(trDt1, type="tip")[rownames(dtTest1), ], dtTest1)
    trDt2 <- readNexus(file=treepluschar, type="all", levels.uniform=FALSE)
    checkIdentical(labels(trDt2), labTr)   # check labels
    checkIdentical(edgeLength(trDt2), eTr) # check edge lengths
    checkIdentical(nodeType(trDt2), nTtr)  # check node types
    checkIdentical(class(trDt2), p4d)      # check class
    checkIdentical(tdata(trDt2, type="tip")[rownames(dtTest2), ], dtTest2)
}
