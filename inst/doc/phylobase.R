## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(
    fig.keep='none', dev='pdf', fig.width=6, fig.height=6,
    latex.options.color="usenames,dvipsnames"
)

## ----randtree1,fig.keep='none',tidy=FALSE--------------------------------
library(ape)
set.seed(1)  ## set random-number seed
rand_tree <- rcoal(10) ## Make a random tree with 10 tips
plot(rand_tree)

## ----convtree,fig.keep='none'--------------------------------------------
library(phylobase)
# convert rand_tree to a phylo4 object
rand_p4_tree <- as(rand_tree, "phylo4")
plot(rand_p4_tree)

## ----geodata,tidy=FALSE--------------------------------------------------
library(phylobase)
data(geospiza_raw)
## what does it contain?
names(geospiza_raw)

## ----convgeodata---------------------------------------------------------
(g1 <- as(geospiza_raw$tree, "phylo4"))

## ----nodelabelgeodata----------------------------------------------------
nodeLabels(g1)

## ------------------------------------------------------------------------
nodeLabels(g1) <- paste("N", nodeId(g1, "internal"), sep="")
head(g1, 5)

## ----sumgeodata----------------------------------------------------------
summary(g1)

## ----tiplabelgeodata-----------------------------------------------------
tipLabels(g1)

## ----modlabelsgeodata----------------------------------------------------
tipLabels(g1) <- tolower(tipLabels(g1))

## ----nodenumbergeodata---------------------------------------------------
nodeId(g1, type='all')

## ----hasbrlengeodata-----------------------------------------------------
hasEdgeLength(g1)

## ----edgeLength-geodata--------------------------------------------------
edgeLength(g1)

## ----edgelabelgeodata----------------------------------------------------
edgeLabels(g1)

## ----edgelabel-assign-geodata--------------------------------------------
edgeLabels(g1)["23-24"] <- "an edge"
edgeLabels(g1)

## ----getEdge-geodata-----------------------------------------------------
getEdge(g1, 24) # default uses descendant node
getEdge(g1, 24, type="ancestor") # edges using ancestor node

## ----getEdge-edgeLength--------------------------------------------------
edgeLength(g1)[getEdge(g1, 24)]
edgeLength(g1)[getEdge(g1, 24, "ancestor")]

## ----rootedgeodata-------------------------------------------------------
isRooted(g1)

## ----rootnodegeodata-----------------------------------------------------
rootNode(g1)

## ----polygeodata---------------------------------------------------------
hasPoly(g1)

## ----ultrametric-geodata-------------------------------------------------
isUltrametric(g1)

## ----nodeDepth-geodata---------------------------------------------------
nodeDepth(g1, 23)
depthTips(g1)

## ----dataprep------------------------------------------------------------
g1 <- as(geospiza_raw$tree, "phylo4")
geodata <- geospiza_raw$data

## ----echo=FALSE, results='hide'------------------------------------------
geodata <- geospiza_raw$data

## ----echo=FALSE, results='hide'------------------------------------------
g2 <- phylo4d(g1, geodata, missing.data="OK", extra.data="OK")

## ----geomerge3, results='hide'-------------------------------------------
g1sub <- prune(g1, "olivacea")
g1B <- phylo4d(g1sub, geodata)

## ----geomergesum---------------------------------------------------------
summary(g2)

## ----geoextract,results='hide'-------------------------------------------
subset(g2, tips.include=c("fuliginosa", "fortis", "magnirostris",
            "conirostris", "scandens"))
subset(g2, node.subtree=21)
subset(g2, mrca=c("scandens", "fortis"))

## ----geodrop, results='hide'---------------------------------------------
subset(g2, tips.exclude=c("fuliginosa", "fortis", "magnirostris",
            "conirostris", "scandens"))
subset(g2, tips.exclude=names(descendants(g2, MRCA(g2, c("difficilis",
               "fortis")))))


## ----getnode-------------------------------------------------------------
data(geospiza)
getNode(geospiza, 10)
getNode(geospiza, "pauper")

## ----getnode2------------------------------------------------------------
getNode(geospiza)
getNode(geospiza, 10:14)
getNode(geospiza, "melanogaster", missing="OK") # no warning
getNode(geospiza, "melanogaster", missing="warn") # warning!

## ----children------------------------------------------------------------
children(geospiza, 16)
ancestor(geospiza, 16)

## ----descendants---------------------------------------------------------
descendants(geospiza, 16)    # by default returns only the tips
descendants(geospiza, "all") # also include the internal nodes
ancestors(geospiza, 20)
ancestors(geospiza, 20, "ALL") # uppercase ALL includes self

## ----siblings------------------------------------------------------------
siblings(geospiza, 20)
siblings(geospiza, 20, include.self=TRUE)

## ----mrca----------------------------------------------------------------
MRCA(geospiza, 1:6)
shortestPath(geospiza, 4, "pauper")

## ----rtree2--------------------------------------------------------------
set.seed(1001)
tree <- as(rcoal(12), "phylo4")

## ----vcvphylo------------------------------------------------------------
vmat <- as(tree, "phylo4vcov")
vmat <- cov2cor(vmat)
library(MASS)
trvec <- mvrnorm(1, mu=rep(0, 12), Sigma=vmat)

## ----plotvcvphylo--------------------------------------------------------
treed <- phylo4d(tree, tip.data=as.data.frame(trvec))
plot(treed)

