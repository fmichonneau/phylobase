## try to read NEXUS files
library(phylobase)
fn <- system.file("nexusfiles/treepluscharV01.nex",package="phylobase")
td<-NexusToPhylo4D(fn)
summary(td)
## would try plotting, but typically don't have enough room
## to plot data
## Error in .local(x, ...) : 
##    No room left to plot data; please try reducing ratio.tree or cex.label.
plot(as(td,"phylo4"))

## try to read a nexus file where the newick string describing the tree is split
## across several lines
multiLine <- system.file("nexusfiles/MultiLineTrees.nex",package="phylobase")
multiLineTrees <-NexusToPhylo4(multiLine)
summary(multiLineTrees)

