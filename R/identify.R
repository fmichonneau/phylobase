setGeneric("identify")
## FIXME: warning on install because
##   "identify.node.R" is read before "phylo4.R"
##   (class definitions)
setMethod("identify","phylo4",
           function(x,n=1,...) {
               plot(x)
               last <- .last_plot.phylo # information return from plot.phylo
               N.tip <- nTips(x)
               N.node <- nNodes(x)
               tips <- labels(x)
               if (!hasNodeLabels(x)) {
                   nodes <- (N.tip+1):(N.tip+N.node)
               } else nodes <- NodeLabels(x)
               labs<-c(rep("",N.tip), nodes)
               click <- identify(last$xx, last$yy, labels=labs, n=n)
               ##    if (click > N.tip) {
               ##    stop("this is a tip, you have to choose a node\n")
               ##}
               names(click) <- labs[click]
               return(click)
           })
