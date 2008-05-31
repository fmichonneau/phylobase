
NexusToPhylo4 <- function(fileToRead,multi=FALSE) {
#This returns a single phylo4 object if the nexus file has one tree; 
#it returns a list of phylo4 objects if the nexus file has more than one or multi==TRUE
#Note that a list of phylo4 objects is not the same as a multiphylo4d object
	trees<-c("Failure");
	filename<-fileToRead;
	params <- list(filename=fileToRead)
	
# Check that params is properly formatted.
	if(!is.list(params) || length(params) == 0) {
		stop("The params parameter must be a non-empty list");
	}
	
# Finally ready to make the call...
    intreesstring <- .Call("ReadTreesWithNCL", params,
						   PACKAGE="phylobase")
	
	intreesphylolist <- read.nexustreestring(intreesstring);
	if (length(intreesphylolist)>1 || multi) {
		trees<-list()
		for (i in 1:length(intreesphylolist)) {
			trees[[i]]<-as(intreesphylolist[[i]], "phylo4");
		}
	}
	else {
		trees<-as(intreesphylolist[[1]], "phylo4");
	}
	trees
}

NexusToDataFrame <- function(fileToRead,allchar=FALSE, polymorphictomissing=TRUE, levelsall=TRUE) {
#This returns a single phylo4D object if the nexus file has one tree; 
#it returns a list of phylo4D objects if the nexus file has more than one or multi==TRUE
#Note that a list of phylo4D objects is not the same as a multiphylo4d object
	output<-c("Failure")
	filename<-fileToRead
	params <- list(filename=fileToRead, allchar=allchar, polymorphictomissing=polymorphictomissing, levelsall=levelsall)
	
# Check that params is properly formatted.
	if(!is.list(params) || length(params) == 0) {
		stop("The params parameter must be a non-empty list")
	}
	
	incharsstring <- .Call("ReadCharsWithNCL",params,
						   PACKAGE="phylobase")
	print(incharsstring)
	tipdata<-eval(parse(text=incharsstring))
	tipdata
}


NexusToPhylo4D <- function(fileToRead,multi=FALSE,allchar=FALSE, polymorphictomissing=TRUE, levelsall=TRUE) {
#This returns a single phylo4D object if the nexus file has one tree; 
#it returns a list of phylo4D objects if the nexus file has more than one or multi==TRUE
#Note that a list of phylo4D objects is not the same as a multiphylo4d object
	output<-c("Failure")
	tipdata<-NexusToDataFrame(fileToRead,allchar, polymorphictomissing, levelsall);
	intreesphylolist<-NexusToPhylo4(fileToRead,multi=TRUE);
	if (length(intreesphylolist)>1 || multi) {
		output<-list()
		for (i in 1:length(intreesphylolist)) {
			output[[i]]<-phylo4d(as(intreesphylolist[[i]], "phylo4"), tip.data = tipdata)
		}
	}
	else {
		output<-phylo4d(as(intreesphylolist[[1]], "phylo4"), tip.data = tipdata)
	}
	output
}


read.nexustreestring <- function(X)
{
#Returns list of phylo objects (not multi.phylo, and always a list, even if there is only one element
#X is a character vector, each element is one line from a treefile
#This is based almost entirely on read.nexus from APE (Emmanuel Paradis).
	
## first remove all the comments
#	print(X)
#	print(mode(X))
#X<-c(X)
#print(X)
#print(mode(X))
## speedier removal of comments pc 13 April 2008
X <- lapply(X, gsub, pattern = "\\[[^\\]]*\\]", replacement = "")

	X<-unlist(strsplit(unlist(X),c("\n")))
#	print(X)
#	print(mode(X))

## Old comment removal code
    ## LEFT <- grep("\\[", X)
    ## RIGHT <- grep("\\]", X)
    ## if (length(LEFT)) {
    ##     for (i in length(LEFT):1) {
    ##         if (LEFT[i] == RIGHT[i]) {
    ##             X[LEFT[i]] <- gsub("\\[.*\\]", "", X[LEFT[i]])
    ##         } else {
    ##             X[LEFT[i]] <- gsub("\\[.*", "", X[LEFT[i]])
    ##             X[RIGHT[i]] <- gsub(".*\\]", "", X[RIGHT[i]])
    ##             if (LEFT[i] < RIGHT[i] - 1) X <- X[-((LEFT[i] + 1):(RIGHT[i] - 1))]
    ##         }
    ##     }
    ## }

    X <- gsub("ENDBLOCK;", "END;", X, ignore.case = TRUE)
    endblock <- grep("END;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- FALSE
    if (length(i2) == 1) if (i2 > i1) translation <- TRUE
    if (translation) {
        end <- semico[semico > i2][1]
        x <- paste(X[i2:end], sep = "", collapse = "")
        x <- gsub("TRANSLATE", "", x, ignore.case = TRUE)
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[x != ""]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
    }
    ## find the first element in the vector of tree strings which starts with "TREE " or "tree " 
    start <- grep("^\\W?[Tt][Rr][Ee][Ee]\\W", X)[1]
    end <- endblock[endblock > i1][1] - 1
#	print(X)
#	print(mode(X))
    tree <- paste(X[start:end], sep = "", collapse = "")
    tree <- gsub(" ", "", tree)
    tree <- unlist(strsplit(tree, "[=;]"))
    tree <- tree[grep("[\\(\\)]", tree)]
    nb.tree <- length(tree)
    STRING <- as.list(tree)
    trees <- list()
    for (i in 1:nb.tree) {
        ## TODO tree.build and clado.build are ape functions
        obj <- if (length(grep(":", STRING[[i]]))) tree.build(STRING[[i]]) else clado.build(STRING[[i]])
        if (translation) {
            for (j in 1:length(obj$tip.label)) {
                ind <- which(obj$tip.label[j] == TRANS[, 1])
                obj$tip.label[j] <- TRANS[ind, 2]
            }
            if (!is.null(obj$node.label)) {
                for (j in 1:length(obj$node.label)) {
                    ind <- which(obj$node.label[j] == TRANS[, 1])
                    obj$node.label[j] <- TRANS[ind, 2]
                }
            }
        }
        trees[[i]] <- obj
## Check here that the root edge is not incorrectly represented
## in the object of class "phylo" by simply checking that there
## is a bifurcation at the root (node "-1")
        if (sum(trees[[i]]$edge[, 1] == "-1") == 1 && dim(trees[[i]]$edge)[1] > 1) {
            warning("The root edge is apparently not correctly represented\nin your tree: this may be due to an extra pair of\nparentheses in your file; the returned object has been\ncorrected but your file may not be in a valid Newick\nformat")
            ind <- which(trees[[i]]$edge[, 1] == "-1")
            trees[[i]]$root.edge <- trees[[i]]$edge.length[ind]
            trees[[i]]$edge.length <- trees[[i]]$edge.length[-ind]
            trees[[i]]$edge <- trees[[i]]$edge[-ind, ]
            for (j in 1:length(trees[[i]]$edge))
			if (as.numeric(trees[[i]]$edge[j]) < 0)
			trees[[i]]$edge[j] <- as.character(as.numeric(trees[[i]]$edge[j]) + 1)
## Check a second time and if there is still a problem...!!!
            if(sum(trees[[i]]$edge[, 1] == "-1") == 1)
			stop("There is apparently two root edges in your file: cannot read tree file")
        }
    }
	trees
}
