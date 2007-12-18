
## REQUIRED for all trees
check_phylo4 <- function(object) {
	## check consistency of edges etc.
	N <- nrow(object@edge)
	if (hasEdgeLength(object) && length(object@edge.length) != N)
		return("edge lengths do not match number of edges")
	if (length(object@tip.label)+object@Nnode-1 != N)
		return("number of tip labels not consistent with number of edges and nodes")
	return(TRUE)
}


check_data <- function(object,
		use.tip.names=TRUE,
		missing.tip.data=c("fail","OK","warn"),
		extra.tip.data=c("fail","OK","warn"),
		default.tip.names=c("warn","OK","fail"),
		use.node.names=FALSE,
		missing.node.data=c("OK","warn","fail"),		
		extra.node.data=c("OK","warn","fail"),												
		default.node.names=c("warn","OK","fail"),...)							 

{

	## tip default: use names, require names, must match exactly
	#use.tip.names=match.arg(use.tip.names)
	missing.tip.data <- match.arg(missing.tip.data)
	extra.tip.data <- match.arg(extra.tip.data)
	default.tip.names <- match.arg(default.tip.names)
	
	## node default: don't use node names, don't require names, do not need to match exactly
	#use.node.names=match.arg(use.node.names)
	missing.node.data <- match.arg(missing.node.data)
	extra.node.data <- match.arg(extra.node.data)
	missing.node.names <- match.arg(default.node.names)
		
	## for each set of data, check for names, missing and extra data and take appropriate actions
	
	## tip data checks
	## if tip.data exist
	if (!all(dim(object@tip.data)==0)) {
		## if we want to use tip.names
		if (use.tip.names) {
		
		#check for default names
		if (all(row.names(object@tip.data) == 1:length(row.names(object@tip.data)))) {
				#no tip.names
				if (default.tip.names == "fail") {
					stop("Tip data have default names and may not match tree tip labels. Consider using the use.tip.names=FALSE option.")
				}
				else if (default.tip.names == "warn") {
					warning("Tip data have default names and may not match tree tip labels. Consider using the use.tip.names=FALSE option.")
				}
			}
			
			# TODO add check for when length is equal but there is no overlap in names
			## check tip names
			## check for missing or extra tip data (relative to tree taxa)
			if (setequal(row.names(object@tip.data), object@tip.label)) {
				##names are perfect match - ok
				return(TRUE)
			}
			else
			{
				#we know the tree taxa and tip.data taxa are not a perfect match
				
				#if tip.data taxa are subset of tree taxa, check missing.tip.data arg and act accordingly
				if (!all(object@tip.label %in% row.names(object@tip.data))) {
					#we know it's not an exact match - we have missing.tip.data - take action
					if (!any(object@tip.label %in% row.names(object@tip.data))) {
						if (missing.tip.data == "fail") {
							stop("Tip data names do not match tree tip labels.")
						}
						else if (missing.tip.data == "warn") {
							warning("Tip data names do not match tree tip labels.")
						}
					
					}
					else
					{
						if (missing.tip.data == "fail") {
							stop("Tip data names are a subset of tree tip labels.")
						}
						else if (missing.tip.data == "warn") {
							warning("Tip data names are a subset of tree tip labels.")
						}
					}
					#else ok
				}
				
				#if tree taxa are subset of tip.data, check extra.tip arg and act accordingly
				if (!all(row.names(object@tip.data) %in% object@tip.label)) {
					#we know it's not an exact match - we have extra.tip.data - take action
					#fail
					if (extra.tip.data == "fail") {
						stop("Tip data names are a superset of tree tip labels.")
					}
					#warn
					else if (extra.tip.data == "warn") {
						warning("Tip data names are a superset of phylo4 labels.")
					}
					#else ok
				}
				
				return(TRUE)
			} 
		}
		else
		{
			#don't use tip names or attempt to sort - but check to make sure dimensions match
			if (!(phylo4::nTips(object)==dim(object@tip.data)[1])) {
				stop("Ignoring tip data names. Number of tip data do not match number of tree tips.")
			}
		}
	}

	## node data checks
	## if node.data exist
	if (!all(dim(object@node.data)==0)) {
		## if we want to use node.names
		if (use.node.names) {
		
		#check for default names
		#could be 1:nTips or nTips+1:nEdges
		if (all(row.names(object@node.data) == 1:length(row.names(object@node.data))) 
					|| all(row.names(object@node.data) == (phylo4::nTips(object)+1):nEdges(object)))
		{
				#no node.names
				if (default.node.names == "fail") {
					stop("Node data have default names and may not match node labels. Consider using the use.node.names=FALSE option.")
				}
				else if (default.node.names == "warn") {
					warning("Node data have default names and may not match node labels. Consider using the use.node.names=FALSE option.")
				}
			}
			
			## check node names
			## check for missing or extra node data (relative to tree taxa)
			if (setequal(row.names(object@node.data), object@node.label)) {
				##names are perfect match - ok
				return(TRUE)
			}
			else
			{
				#we know the tree taxa and node.data taxa are not a perfect match
				
				#if node.data taxa are subset of tree taxa, check missing.node.data arg and act accordingly
				if (!all(object@node.label %in% row.names(object@node.data))) {
					#we know it's not an exact match - we have missing.node.data - take action
					if (!any(object@node.label %in% row.names(object@node.data))) {
						if (missing.node.data == "fail") {
							stop("Node data names do not match tree node labels.")
						}
						else if (missing.node.data == "warn") {
							warning("Node data names do not match tree node labels.")
						}
					
					}
					else
					{
						if (missing.node.data == "fail") {
							stop("Node data names are a subset of tree node labels.")
						}
						else if (missing.node.data == "warn") {
							warning("Node data are a subset of phylo4 node labels.")
						}
					}
					#else ok
				}
				
				#if tree taxa are subset of node.data, check extra.node arg and act accordingly
				if (!all(row.names(object@node.data) %in% object@node.label)) {
					#we know it's not an exact match - we have missing.node.data - take action
					#fail
					if (extra.node.data == "fail") {
						stop("Node data names are a superset of tree node labels.")
					}
					#warn
					else if (extra.node.data == "warn") {
						warning("Node data names are a superset of tree node labels.")
					}
					#else ok
				}
				
				return(TRUE)
			} 
		}
		else
		{
			#don't use node names or attempt to sort - but check to make sure dimensions match
			if (!(phylo4::nNodes(object)==dim(object@node.data)[1])) {
				stop("Ignoring node data names. Number of node data do not match number of tree nodes.")
			}
		}
	}

}

attach_data <- function(object,
		use.tip.names=TRUE,
		use.node.names=FALSE,
		...)							 
{
	
	## assumes data have already been checked by check_data!
	
	## clean empty names - if data are imported as 'all' but only tips or nodes have names, clean up
	#if (all(row.names(object@tip.data)==""))
	#	row.names(object@tip.data) <- NULL
	#if (all(row.names(object@node.data)==""))
	#	row.names(object@node.data) <- NULL
	
	## for each set of data, take appropriate actions
	
	## tip data operations:
	## if tip.data exist
	if (!all(dim(object@tip.data)==0)) {
		## if we want to use tip.names
		if (use.tip.names) {
			object@tip.data <- object@tip.data[match(object@tip.label,row.names(object@tip.data)),,drop=FALSE]
		}
		row.names(object@tip.data) <- object@tip.label
	}
	
	## node data operations
	if (!all(dim(object@node.data)==0)) {
		## if we want to use tip.names
		if (use.node.names) {
			object@node.data <- object@node.data[match(object@node.label,row.names(object@node.data)),,drop=FALSE]
		}
		row.names(object@node.data) <- object@node.label
	}
	
	return(object)
	
}
