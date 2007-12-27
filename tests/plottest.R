library(phylobase)
library(geiger)
library(ape)

data(geospiza)
g1 <- as(geospiza$geospiza.tree,"phylo4")
g2 <- phylo4d(g1,geospiza$geospiza.data,missing.tip.data="OK")

par(mfrow=c(1,2))
plot(g1,show.node.label=TRUE)
plot(g2,show.node.label=TRUE)

##  Note the numbering differences!

 g3 = subset(g2,tips.exclude=c("fuliginosa","fortis","magnirostris",
            "conirostris","scandens"))
 plot(extract.tree(g3))  ## phylo4
t1 <- try(plot(g3))                 ## phylo4d -- error
