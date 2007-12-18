################
# subset phylo4
################
#
setMethod("subset", "phylo4",
  function(x,tips.include=NULL,tips.exclude=NULL,mrca=NULL,node.subtree=NULL,...) {

    if (!is.null(tips.include)) {
      if (is.numeric(tips.include)) {
        tips.include <- x@tip.label[tips.include]
      }
      return(prune(x,x@tip.label[!(x@tip.label %in% tips.include)]))
    }
    
    if (!is.null(tips.exclude)) {
      return(prune(x,tips.exclude))
    }
    
    if (!is.null(node.subtree)) {
      return("Method not impelemented yet.")
    }
    
    if (!is.null(mrca)) {
      
      if (length(mrca) != 2) {
        return(paste("MRCA subset requires 2 terminal taxa, you supplied",length(mrca)))
      }
      
      if (is.numeric(mrca)) {
        mrca <- x@tip.label[mrca]
      }
      
     list.nodes<-list(2)
     m<-1:length(x@tip.label)
     
     for (i in 1:2) {
        fils <- m[x@tip.label[]==mrca[i]]
        
        res<-NULL
        repeat
          {
          pere<-x$edge[,1][x$edge[,2]==fils]
          res<-c(res,pere)
          fils<-pere
          if(pere==(length(x@tip.label)+1)) break
          }
        list.nodes[[i]]<-res
       } 
       
       MRCA<-max(intersect(list.nodes[[1]],list.nodes[[2]]))     
       
       fils<-NULL
       pere<-res <- MRCA
       if (MRCA==length(x@tip.label)+1) return(x)
       else {
            repeat
	       {
	    for (i in 1: length(pere)) fils<-c(fils, x@edge[,2][x@edge[,1]==pere[i]])
	    res<-c(res, fils)
        pere<-fils
	    fils<-NULL
	    if (length(pere)==0) break
	       }
       
       tips.exclude <- setdiff(x@tip.label,x@tip.label[res[res<length(x@tip.label)+1]])
       print(tips.exclude)
      return(prune(x,tips.exclude))
      } 
         }
    
    return(x)
  })



setMethod("subset", "phylo", function(x,tips.include=NULL,tips.exclude=NULL,mrca=NULL,node.subtree=NULL,...) {

  x <- as(x,"phylo")
  res <- subset(x,tips.include=NULL,tips.exclude=NULL,mrca=NULL,node.subtree=NULL,...)
  return(res)
  })



#
# PLEASE CODE ME 
#
setMethod("subset", "phylo4d", function(x,tips.include=NULL,tips.exclude=NULL,mrca=NULL,node.subtree=NULL,...) {

# place commercial here.
})



