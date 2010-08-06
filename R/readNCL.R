### This file contains the source code for the functions:
###  - readNCL (generic function)
###  - readNexus (wrapper for readNCL importing Nexus files)
###  - readNewick (wrapper for readNCL importing Newick files)

readNCL <- function(file, simplify=FALSE, type=c("all", "tree", "data"),
                    char.all=FALSE, polymorphic.convert=TRUE,
                    levels.uniform=FALSE, quiet=TRUE,
                    check.node.labels=c("keep", "drop", "asdata"),
                    return.labels=TRUE, file.format=c("nexus", "newick"), ...) {


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
 ncl <- .Call("GetNCL", fileName, parameters, PACKAGE="phylobase")

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
   tipData <- data.frame(tipData)
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
   for (i in 1:length(ncl$trees)) {
     if (length(grep(":", ncl$trees[i]))) {
       listTrees[[i]] <- tree.build(ncl$trees[i])
     }
     else {
       listTrees[[i]] <- clado.build(ncl$trees[i])
     }
   }
   listTrees <- lapply(listTrees, function(tr) {       
     if (length(ncl$taxaNames) == nTips(tr)) {
       tr$tip.label <- ncl$taxaNames[as.numeric(tr$tip.label)]     
     }
     else stop("phylobase doesn't deal with multiple taxa block at this time.")
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
                       return.labels=TRUE, ...) {

  return(readNCL(file=file, simplify=simplify, type=type, char.all=char.all,
          polymorphic.convert=polymorphic.convert, levels.uniform=levels.uniform,
          quiet=quiet, check.node.labels=check.node.labels,
          return.labels=return.labels, file.format="nexus", ...))
}

readNewick <- function(file, simplify=FALSE, quiet=TRUE,
                       check.node.labels=c("keep", "drop", "asdata"), ...) {

  return(readNCL(file=file, simplify=simplify, quiet=quiet,
                 check.node.labels=check.node.labels, file.format="newick", ...))
}
