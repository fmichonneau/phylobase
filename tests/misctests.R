library(phylobase)
library(ape)

data(geospiza)
geospiza0 <-
  list(geospiza.tree=as(geospiza,"phylo"),geospiza.data=tdata(geospiza))
## push data back into list form as in geiger

t1 <-  try(p1 <- phylo4d(geospiza0$geospiza.tree,geospiza0$geospiza.data))
## Error in check_data(res, ...) : 
##   Tip data names are a subset of tree tip labels.

p2 <- as(geospiza0$geospiza.tree,"phylo4")
plot(p2)

lab1 <- labels(p2)
lab2 <- rownames(geospiza0$geospiza.data)

lab1[!lab1 %in% lab2]  ## missing data
lab2[!lab2 %in% lab1]  ## extra data (none)
p1 <- phylo4d(p2,geospiza0$geospiza.data,missing.tip.data="warn")
p1 <- phylo4d(p2,geospiza0$geospiza.data,missing.tip.data="OK")

plot(p1)
plot(p1,show.node.label=TRUE)
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
q1 <- phylo4d(r1,tip.data=tipdat,node.data=data.frame(a=6:9))
q2 <- prune(q1,1)
summary(q2)

tipdat2 <- tipdat
row.names(tipdat2)[1] <- "s1"
t1 <- try(q1 <- phylo4d(r1,tip.data=tipdat2))

plot(q2)
plot(q2,type="cladogram")
## plot(p2,type="dotchart",labels.nodes=nodeLabels(p2))
## trace("plot", browser, signature = c("phylo4d","missing"))
labels(q1) <- paste("q",1:5,sep="")
nodeLabels(q1) <- paste("n",1:4,sep="")
p3 <- phylo4d(r1,tip.data=tipdat,node.data=data.frame(b=6:9))
summary(p3)

plot(p1)

plot(subset(p1,tips.include=c("fuliginosa","fortis","magnirostris",
            "conirostris","scandens")))
## better error?
## Error in phy$edge[, 2] : incorrect number of dimensions

if(dev.cur() == 1) get(getOption("device"))()
plot(subset(p2,tips.include=c("fuliginosa","fortis","magnirostris",
            "conirostris","scandens")))

plot(p2,show.node.label=TRUE)

library(ape)
example(read.tree)

z <- as(tree.owls,"phylo4")

example("phylo4d")
obj1 <- obj2 <- obj3 <- phylo4d(as(tree.owls,"phylo4"),data.frame(wing=1:4,color=factor(c("b","w","b","b")), tail=runif(4)*10), use.tip.names=FALSE)

obj2@tip.data <- as.data.frame(obj2@tip.data[,1])
obj3@tip.data <- cbind(obj1@tip.data,obj2@tip.data)
obj4 <- obj1
obj4$tip.data[2,3] <- NA
obj4$tip.data[1,1] <- NA

obj4@node.label <- character(0)

obj5 <- obj1
tdata(obj4) <- subset(tdata(obj4),select=sapply(tdata(obj4),class)=="numeric")

treePlot(obj4)
