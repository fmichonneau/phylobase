library(phylobase)
library(ape)
r1 <- rcoal(5)

## trace("phylo4d", browser, signature = "phylo")
## untrace("phylo4d", signature = "phylo")
tipdat <- data.frame(a=1:5,row.names=r1$tip.label)
p1 <- phylo4d(r1,tip.data=tipdat,node.data=data.frame(a=6:9))
p2 <- prune(p1,1)
summary(p2)

labels(p1) <- paste("q",1:5,sep="")
NodeLabels(p1) <- paste("n",1:4,sep="")
p3 <- phylo4d(r1,tip.data=tipdat,node.data=data.frame(b=6:9))
summary(p3)

plot(p1)
