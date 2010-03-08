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
  # function (file, simplify=TRUE, type=c("all", "tree", "data"),
  #   char.all=FALSE, polymorphic.convert=TRUE, levels.uniform=TRUE,
  #   check.node.labels=c("keep", "drop", "asdata"))
 co1 <- readNexus(file=co1File, check.node.labels="asdata")
 co1Tree1 <- co1[[1]]
 co1Tree2 <- co1[[2]]
 ## Check labels
 labCo1 <- c("Cow", "Seal", "Carp", "Loach", "Frog", "Chicken", "Human",
             "Mouse", "Rat", "Whale", NA, NA, NA, NA, NA, NA, NA, NA)
 names(labCo1) <- 1:18
 checkIdentical(labels(co1Tree1), labCo1)
 checkIdentical(labels(co1Tree2), labCo1)
 ## Check edge lengths
 eLco1 <- c(0.143336, 0.225087, 0.047441, 0.055934, 0.124549, 0.204809, 0.073060, 0.194575,
            0.171296, 0.222039, 0.237101, 0.546258, 0.533183, 0.154442, 0.134574, 0.113163,
            0.145592)
 names(eLco1) <- c("11-1", "11-2", "11-12", "12-13", "13-14", "14-15", "15-16", "16-17", "17-3",
                   "17-4", "16-5", "15-6", "14-7", "13-18", "18-8", "18-9", "12-10")
 checkIdentical(edgeLength(co1Tree1), eLco1)
 checkIdentical(edgeLength(co1Tree2), eLco1)
 ## Check node type
 nTco1 <-  c("tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip",
             "tip", "internal", "internal", "internal", "internal", "internal",
             "internal", "internal", "internal")
 names(nTco1) <- 1:18
 checkIdentical(nodeType(co1Tree1), nTco1)
 checkIdentical(nodeType(co1Tree2), nTco1)
 ## Check label values
 lVco1 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.93, 0.88, 0.99, 1.00,
            0.76, 1.00, 1.00)
 checkIdentical(as(co1Tree1, "data.frame")$labelValues, lVco1)

 ## Mutli Lines
 multiLines <- readNexus(file=multiLinesFile)
 ## load correct representation and make sure that the trees read
 ## match it
 load(mlFile)
 checkIdentical(multiLines[[1]], ml1)
 checkIdentical(multiLines[[2]], ml1)
 rm(ml)
}
