library(phylobase)

set.seed(1)
t0 <- rcoal(5)
t0$edge

plot(t0)

t1<-as(t0,"phylo4")
t5 <- as(t1,"phylo")
stopifnot(identical(t0,t5))

## t2<-as(t1,"phylo4vcov")
## t3<-as(t2,"phylo4")
## t4<-as(t3,"phylo")

## plot(test.tree) #CRASHES R!
