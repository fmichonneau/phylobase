library(phylobase)
library(ape)
library(geiger)

data(geospiza)
t1 <-  try(p1 <- phylo4d(geospiza$geospiza.tree,geospiza$geospiza.data))
## Error in check_data(res, ...) : 
##   Tip data names are a subset of tree tip labels.

p2 <- as(geospiza$geospiza.tree,"phylo4")
plot(p2)

lab1 <- labels(p2)
lab2 <- rownames(geospiza$geospiza.data)

lab1[!lab1 %in% lab2]  ## missing data
lab2[!lab2 %in% lab1]  ## extra data (none)
p1 <- phylo4d(p2,geospiza$geospiza.data,missing.tip.data="OK")
## plot(p1) ## fails -- because of missing data?

## one way to deal with it:

p1B <- prune(p1,tip="olivacea")

## or ...
p1C <- na.omit(p1)

labels(p1C) <- tolower(labels(p1C))

## trace("prune",browser,signature="phylo4d")
r1 <- rcoal(5)

## trace("phylo4d", browser, signature = "phylo")
## untrace("phylo4d", signature = "phylo")
tipdat <- data.frame(a=1:5,row.names=r1$tip.label)
p1 <- phylo4d(r1,tip.data=tipdat,node.data=data.frame(a=6:9))
p2 <- prune(p1,1)
summary(p2)

plot(p2)
plot(p2,type="dotchart")
## plot(p2,type="dotchart",labels.nodes=NodeLabels(p2))
## trace("plot", browser, signature = c("phylo4d","missing"))
labels(p1) <- paste("q",1:5,sep="")
NodeLabels(p1) <- paste("n",1:4,sep="")
p3 <- phylo4d(r1,tip.data=tipdat,node.data=data.frame(b=6:9))
summary(p3)

plot(p1)

