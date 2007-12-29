library(phylobase)
library(geiger)
library(ape)

data(geospiza)
g1 <- as(geospiza$geospiza.tree,"phylo4")
g2 <- phylo4d(g1,geospiza$geospiza.data,missing.tip.data="OK")

par(mfrow=c(1,2))
plot(g1,show.node.label=TRUE)
plot(g2,show.node.label=TRUE)

g2B <- as(g2,"phylog")
##  Note the numbering differences!

## round trip 
g2C <- as(read.tree(text=write.tree(as(g1,"phylo"))),"phylo4")
## comes back in same order
plot(g1,show.node.label=TRUE)
plot(g2C,show.node.label=TRUE)

g3 = subset(g2,tips.exclude=c("fuliginosa","fortis","magnirostris",
                 "conirostris","scandens"))
plot(extract.tree(g3))  ## phylo4
t1 <- try(plot(g3))                 ## phylo4d -- error

## Playing with new ways of plotting

dist1 <- cophenetic.phylo(as(g2,"phylo"))
mdspos <- isoMDS(dist1)$points
par(mfrow=c(2,2))
plot(g1)
plot(mdspos,type="n")
text(mdspos[,1],mdspos[,2],abbreviate(rownames(mdspos)))
cmdpos <- cmdscale(dist1)
plot(cmdpos,type="n")
text(cmdpos[,1],cmdpos[,2],abbreviate(rownames(mdspos)))

## never mind, I don't know how to construct a useful
##  2D color space anyway ...
