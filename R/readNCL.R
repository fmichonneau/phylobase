### This file contains the source code for the functions:
###  - readNCL (generic function)
###  - readNexus (wrapper for readNCL importing Nexus files)
###  - readNewick (wrapper for readNCL importing Newick files)

#' Create a phylo4, phylo4d or data.frame object from a Nexus or a Newick file
#' 
#' \code{readNexus} reads a Nexus file and outputs a \code{phylo4} or
#' \code{phylo4d} or \code{data.frame} object.
#' 
#' \code{readNewick} reads a Newick file and outputs a \code{phylo4} or
#' \code{phylo4d} object.
#' 
#' 
#' \code{readNexus} extracts data held in a Nexus file, specifically from DATA,
#' CHARACTER or TREES blocks present in the file. The \code{type} argument
#' specifies which of these is returned: \describe{ \item{data}{will only
#' return a \code{data.frame} of the contents of all DATA and CHARACTER
#' blocks.} \item{tree}{will only return a \code{phylo4} object of the contents
#' of the TREES block.} \item{all}{if only data or a tree are present in the
#' file, this option will act as the options above, returning either a
#' \code{data.frame} or a \code{phylo4} object respectively. If both are
#' present then a \code{phylo4d} object is returned containing both.} } The
#' function returns \code{NULL} if the \code{type} of data requested is not
#' present in the file, or if neither data nor tree blocks are present.
#' 
#' Depending on the context \code{readNexus} will call either the \code{phylo4}
#' or \code{phylo4d} constructor. In addition with \code{type="all"}, the
#' \code{phylo4d} constructor will be used if
#' \code{check.node.labels="asdata"}.
#' 
#' \code{readNewick} imports newick formatted tree files and will return a
#' \code{phylo4} or a \code{phylo4d} object if the option
#' \code{check.node.labels="asdata"} is invoked.
#' 
#' For both \code{readNexus} and \code{readNewick}, the options for
#' \code{check.node.labels} can take the values: \describe{ \item{keep}{the
#' node labels of the trees will be passed as node labels in the \code{phylo4}
#' object} \item{drop}{the node labels of the trees will be ignored in the
#' \code{phylo4} object} \item{asdata}{the node labels will be passed as data
#' and a \code{phylo4d} object will be returned.} } If you use the option
#' \code{asdata} on a file with no node labels, a warning message is issued,
#' and thus \code{check.node.labels} takes the value \code{drop}.
#' 
#' For both \code{readNexus} and \code{readNewick}, additional arguments can be
#' passed to the constructors such as \code{annote}, \code{missing.data} or
#' \code{extra.data}. See the documentation of \link{phylo4-methods},
#' \link{phylo4d} and \link{formatData} for the complete list of options.
#' 
#' @name Import Nexus and Newick files
#' @aliases readNexus readNCL readNewick
#' @docType methods
#' @param file a Nexus file for \code{readNexus} or a file that contains Newick
#' formatted trees for \code{readNewick}
#' @param simplify If there are multiple trees in the file, only the first one
#' is returned if TRUE and a list of phylo4/phylo4d objects is returned if the
#' file contains multiple trees.
#' @param type Determines which type of objects to return, if present in the
#' file (see Details).
#' @param char.all If TRUE, returns all characters, even those excluded in the
#' NEXUS file
#' @param polymorphic.convert If TRUE, converts polymorphic characters to
#' missing data
#' @param levels.uniform If TRUE, uses the same levels for all characters
#' @param quiet If FALSE the output of the NCL interface is printed. This is
#' mainly for debugging purposes. This option can considerably slow down the
#' process if the tree is big or there are many trees in the file.
#' @param check.node.labels Determines how the node labels in the Nexus or
#' Newick files should be treated in the phylo4 object, see Details for more
#' information.
#' @param return.labels Determines whether state names (if TRUE) or state codes
#' should be returned.
#' @param check.names logical. If \sQuote{TRUE} then the names of the
#' characters from the NEXUS file are checked to ensure that they are
#' syntactically valid variable names and are not duplicated.  If necessary
#' they are adjusted (by \sQuote{make.names}) so that they are.
#' @param convert.edge.length logical. If \sQuote{TRUE} negative edge lengths
#' are replaced with 0. At this time \code{phylobase} does not accept objects
#' with negative branch lengths, this workaround allows to still use trees with
#' negative branch lengths are an artifact of the method used to build the
#' tree.
#' @param \dots Additional arguments to be passed to phylo4 or phylo4d
#' constructor (see Details)
#' @return Depending on the value of \code{type} and the contents of the file,
#' one of: a \code{data.frame}, a \linkS4class{phylo4} object, a
#' \linkS4class{phylo4d} object or \code{NULL}.  If several trees are included
#' in the Nexus file and the option \code{simplify=FALSE} a list of
#' \linkS4class{phylo4} or \linkS4class{phylo4d} objects is returned.
#' @note Underscores in state labels (i.e. trait or taxon names) will be
#' translated to spaces when read by NCL. Unless \code{check.names=FALSE},
#' trait names will be converted to valid R names (see
#' \code{\link{make.names}}) on input to R, so spaces will be translated to
#' periods.
#' @author Brian O'Meara, Francois Michonneau, Derrick Zwickl
#' @seealso the \linkS4class{phylo4d} class, the \linkS4class{phylo4} class
#' @keywords misc
readNCL <- function(file, simplify=FALSE, type=c("all", "tree", "data"),
                    char.all=FALSE, polymorphic.convert=TRUE,
                    levels.uniform=FALSE, quiet=TRUE,
                    check.node.labels=c("keep", "drop", "asdata"),
                    return.labels=TRUE, file.format=c("nexus", "newick"),
                    check.names=TRUE, convert.edge.length=FALSE, ...) {

  ## turn on to TRUE to test new way of building trees in NCL
 experimental <- FALSE

 file <- path.expand(file)
 type <- match.arg(type)
 check.node.labels <- match.arg(check.node.labels)
 file.format <- match.arg(file.format)
 if (file.format == "newick") file.format <- "relaxedphyliptree" 
 
 if (type == "all" || type == "data") {
   returnData <- TRUE
 }
 else {
   returnData <- FALSE
 }
 if (type == "all" || type == "tree") {
   returnTrees <- TRUE
 }
 else {
   returnTrees <- FALSE
 }

 fileName <- list(fileName=file, fileFormat=file.format)
 parameters <- c(char.all, polymorphic.convert, levels.uniform, returnTrees, returnData)

 ## GetNCL returns a list containing:
 ##  $taxaNames: names of the taxa (from taxa block, implied or declared)
 ##  $treeNames: the names of the trees
 ##  $trees: a vector of (untranslated) Newick strings
 ##  $dataTypes: data type for each character block of the nexus file (length = number of chr blocks)
 ##  $nbCharacters: number of characters in each block (length = number of chr blocks)
 ##  $charLabels: the labels for the characters, i.e. the headers of the data frame to be returned
 ##    (length = number of chr blocks * sum of number of characters in each block)
 ##  $nbStates: the number of states of each character (equals 0 for non-standard types, length = number
 ##    of characters)
 ##  $stateLabels: the labels for the states of the characters, i.e. the levels of the factors to be returned
 ##  $dataChr: string that contains the data to be returned
 ncl <- GetNCL(fileName, parameters)

 ## Return Error message
 if (length(ncl) == 1 && names(ncl) == "ErrorMsg") {
   stop(ncl$ErrorMsg)
 }
 
 if (!quiet) print(ncl)

 ## Disclaimer
 if (!length(grep("\\{", ncl$dataChr)) && return.labels && !polymorphic.convert) {
   stop("At this stage, it's not possible to use the combination: ",
        "return.labels=TRUE and polymorphic.convert=FALSE for datasets ",
        "that contain polymorphic characters.")
 }
 
 if (returnData && length(ncl$dataChr)) {
   tipData <- vector("list", length(ncl$dataChr))
   for (iBlock in 1:length(ncl$dataTypes)) {
     chrCounter <- ifelse(iBlock == 1, 0, sum(ncl$nbCharacters[1:(iBlock-1)]))
     if (ncl$dataTypes[iBlock] == "Continuous") {       
       for (iChar in 1:ncl$nbCharacters[iBlock]) {
         i <- chrCounter + iChar
         tipData[[i]] <- eval(parse(text=ncl$dataChr[i]))
         names(tipData)[i] <- ncl$charLabels[i]
       }
     }
     else {
 
       if (ncl$dataTypes[iBlock] == "Standard") {
         iForBlock <- integer(0)
         for (iChar in 1:ncl$nbCharacters[iBlock]) {      
           i <- chrCounter + iChar
           iForBlock <- c(iForBlock, i)
           lblCounterMin <- ifelse(i == 1, 1, sum(ncl$nbStates[1:(i-1)]) + 1)
           lblCounter <- seq(lblCounterMin, length.out=ncl$nbStates[i])
           tipData[[i]] <- eval(parse(text=ncl$dataChr[i]))
           names(tipData)[i] <- ncl$charLabels[i]
           tipData[[i]] <- as.factor(tipData[[i]])
          
           lbl <- ncl$stateLabels[lblCounter]
           if (return.labels) {
             if (any(nchar(gsub(" ", "", lbl)) == 0)) {
               warning("state labels are missing for \'", ncl$charLabels[i],
                       "\', the option return.labels is thus ignored.")
             }
             else {
               levels(tipData[[i]]) <- lbl
             }
           }          
         }
         if (levels.uniform) {
           allLevels <- character(0)
           for (j in iForBlock) {
             allLevels <- union(allLevels, levels(tipData[[j]]))
           }
           for (j in iForBlock) {
             levels(tipData[[j]]) <- allLevels
           }
         }
       }
       else {
         warning("This datatype is not currently supported by phylobase")
         next
         ## FIXME: different datatypes in a same file isn't going to work
       }
     }
   }
   tipData <- data.frame(tipData, check.names=check.names)
   if (length(ncl$taxaNames) == nrow(tipData)) {
     rownames(tipData) <- ncl$taxaNames
   }
   else stop("phylobase doesn't deal with multiple taxa block at this time.")
 }
 else {
   tipData <- NULL
 }

 if (returnTrees && length(ncl$trees) > 0) {
   listTrees <- vector("list", length(ncl$trees))
   
   if (!experimental) {
     for (i in 1:length(ncl$trees)) {
       if (length(grep(":", ncl$trees[i]))) {
         ## remove comments from the string, they are not being removed if they are part of the node labels
         tStr <- gsub("\\[.*?\\]", "", ncl$trees[i])
         listTrees[[i]] <- tree.build(tStr)
       }
       else {
         tStr <- gsub("\\[.*?\\]", "", ncl$trees[i])
         listTrees[[i]] <- clado.build(tStr)
       }
     }
     listTrees <- lapply(listTrees, function(tr) {       
       if (length(ncl$taxaNames) == nTips(tr)) {
         tr$tip.label <- ncl$taxaNames[as.numeric(tr$tip.label)]     
       }
       else stop("phylobase doesn't deal with multiple taxa block at this time.")
       if (convert.edge.length) {
         tr$edge.length[tr$edge.length < 0] <- 0
       }
       if (is.null(tr$node.label)) {
         if (check.node.labels == "asdata") {
           warning("Could not use value \"asdata\" for ",
                   "check.node.labels because there are no ",
                   "labels associated with the tree")
           check.node.labels <- "drop"
         }
         tr <- phylo4(tr, check.node.labels=check.node.labels, ...)       
       }
       else {
         if (check.node.labels == "asdata") {
           tr <- phylo4d(tr, check.node.labels=check.node.labels, ...)
         }
         else {
           tr <- phylo4(tr, check.node.labels=check.node.labels, ...)
         }
       }
     })
     if (length(listTrees) == 1 || simplify)
       listTrees <- listTrees[[1]]
   }
   else {
     edgeMat <- cbind(ncl$parentVector, c(1:length(ncl$parentVector)))
     edgeLgth <- ncl$branchLengthVector
     edgeLgth[edgeLgth == -1] <- NA
     if (convert.edge.length) {
       edgeLgth[edgeLgth < 0] <- 0
     }
     if (length(ncl$taxaNames) != min(ncl$parentVector)-1) {
       stop("phylobase doesn't deal with multiple taxa block at this time.")
     }
     ## TODO: code node labels in GetNCL
     tr <- phylo4(x=edgeMat, edge.length=edgeLgth, tip.label=ncl$taxaNames,
                  ...)
   }
 }
 else {
   listTrees <- NULL
 }
 
 ###
 switch(type,
        "data" = {
          if (is.null(tipData)) {
            toRet <- NULL
          }
          else {
            toRet <- tipData
          }
        },
        "tree" = {
          if (is.null(listTrees)) {
            toRet <- NULL
          }
          else {           
            toRet <- listTrees                                          
          }
        },
        "all" = {
          if (is.null(tipData) && is.null(listTrees)) {
            toRet <- NULL
          }
          else if (is.null(tipData)) {              
            toRet <- listTrees
          }
          else if (is.null(listTrees)) {
            toRet <- tipData
          }
          else {
            if (length(listTrees) > 1) {              
              toRet <- lapply(listTrees, function(tr)
                              addData(tr, tip.data=tipData, ...))
            }
            else toRet <- addData(listTrees, tip.data=tipData, ...)
          }
        })
 toRet
}
 
readNexus <- function (file, simplify=FALSE, type=c("all", "tree", "data"),
                       char.all=FALSE, polymorphic.convert=TRUE,
                       levels.uniform=FALSE, quiet=TRUE,
                       check.node.labels=c("keep", "drop", "asdata"),
                       return.labels=TRUE, check.names=TRUE, convert.edge.length=FALSE,
                       ...) {

  return(readNCL(file=file, simplify=simplify, type=type, char.all=char.all,
          polymorphic.convert=polymorphic.convert, levels.uniform=levels.uniform,
          quiet=quiet, check.node.labels=check.node.labels,
          return.labels=return.labels, file.format="nexus",
          check.names=check.names, convert.edge.length=convert.edge.length, ...))
}

readNewick <- function(file, simplify=FALSE, quiet=TRUE,
                       check.node.labels=c("keep", "drop", "asdata"),
                       convert.edge.length=FALSE, ...) {

  return(readNCL(file=file, simplify=simplify, quiet=quiet,
                 check.node.labels=check.node.labels, file.format="newick",
                 convert.edge.length=convert.edge.length, ...))
}
