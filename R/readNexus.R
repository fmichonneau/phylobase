readNexus <- function (file, simplify=FALSE, type=c("all", "tree", "data"),
                       char.all=FALSE, polymorphic.convert=TRUE,
                       levels.uniform=TRUE, quiet=TRUE,
                       check.node.labels=c("keep", "drop", "asdata"),
                       return.labels=TRUE, ...) {

    ## file = input nexus file
    ## simplify = if TRUE only keeps the first tree, if several trees are found in
    ##            the Nexus file
    ## type = specify whether to return trees+data as phylo4d object ("all") if
    ##        both are found, returning a data.frame or phylo4 object if only one
    ##        is found, "tree": return a phylo4 object only, regardless of
    ##        whether there are data, "data": return a data.frame (no tree), even
    ##        if a tree is present
    ## char.all = if TRUE, includes even excluded chars in the nexus file
    ## polymorphic.convert = if TRUE, convert polymorphic characters to missing
    ##                       characters
    ## levels.uniform = if TRUE, categorical data are loaded with the same levels,
    ##                  even if one character is missing a state
    ## return.labels = if TRUE, returns the names of the states instead of the
    ##                 the internal codes
    ## quiet = if TRUE, returns the object without printing tree strings (printing
    ##         makes readNexus very slow in the cases of very big trees)
    ## check.node.labels = how to deal with node labels, to be passed to phylo4d
    ##                     constructor

    type <- match.arg(type)
    check.node.labels <- match.arg(check.node.labels)

    output <- c("Failure")
    if (type == "all" || type == "data") {
        params <- list(filename=file, allchar=char.all,
                       polymorphictomissing=polymorphic.convert,
                       levelsall=levels.uniform,
                       returnlabels=return.labels)

        ## Check that params is properly formatted.
        if(!is.list(params) || length(params) == 0) {
            stop("The params parameter must be a non-empty list")
        }
        incharsstring <- .Call("ReadCharsWithNCL", params,
                               PACKAGE="phylobase")
        if (length(incharsstring) > 0) {
            incharsstring <- unlist(strsplit(incharsstring$charstring, "\\|"))
            incharsstring <- incharsstring[nzchar(incharsstring)]

            if (!quiet) print(incharsstring)   # display character string if quiet is FALSE

            iDtType <- seq(from=1, to=length(incharsstring), by=2)
            iCharStrg <- seq(from=2, to=length(incharsstring), by=2)

            datatype <- incharsstring[iDtType]
            charString <- incharsstring[iCharStrg]

            tipdata <- list()
            for (i in 1:length(charString)) {
                if (datatype[i] == "Standard") {
                    ## Remove empty labels for factors
                    charString[i] <- gsub("\\\"\\\"", "", charString[i])
                    charString[i] <- gsub(",+)", ")", charString[i])

                    ## For now, we can't deal with polymorphic characters and their labels
                    if (length(grep("\\{", charString[i])) > 0 &&
                        return.labels) {
                        stop("At this stage, it's not possible to use the combination: ",
                             "return.labels=TRUE for datasets that contain polymorphic ",
                             "characters.")
                    }

                    ## Convert the string to data frame
                    tipdata[[i]] <- eval(parse(text=charString[i]))

                    ## if levels.uniform=TRUE apply the same levels to all characters
                    if (levels.uniform && length(tipdata[[i]]) > 0) {
                        allLevels <- character(0)
                        for (j in 1:ncol(tipdata[[i]])) {
                            allLevels <- union(allLevels, levels(tipdata[[i]][,j]))
                        }
                        for (j in 1:ncol(tipdata[[i]])) {
                            levels(tipdata[[i]][,j]) <- allLevels
                        }
                    }
                }
                else {
                    ## Just convert string to data frame for other datatype
                    tipdata[[i]] <- eval(parse(text=charString[i]))
                }
            }
            finalTipdata <- tipdata[[1]]
            if (length(tipdata) > 1) {
                for(td in tipdata[-1]) {
                    finalTipdata <- cbind(finalTipdata, td)
                }
            }
            tipdata <- finalTipdata
        }
        else {
            tipdata <- NULL
        }
    }
    if (type == "all" || type == "tree") {
        trees <- c("Failure");
        params <- list(filename=file)

        ## Check that params is properly formatted.
        if(!is.list(params) || length(params) == 0) {
            stop("The params parameter must be a non-empty list");
        }

        ## Finally ready to make the call...
        intreesstring <- .Call("ReadTreesWithNCL", params,
                               PACKAGE="phylobase")
        ## Display the string returned by NCL if quiet=FALSE
        if(!quiet) print(intreesstring)
        intreesphylolist <- read.nexustreestring(intreesstring)
        if (length(intreesphylolist)>1 && !simplify) {
            trees <- list()
            for (i in 1:length(intreesphylolist)) {
                if(identical(check.node.labels, "asdata")) {
                    if(is.null(intreesphylolist[[i]]$node.label)) {
                        warning("Could not use value \"asdata\" for ",
                                "check.node.labels because there are no ",
                                "labels associated with the tree ", i)
                        check.node.labels <- "drop"
                    }
                    trees[[i]] <- phylo4d(intreesphylolist[[i]],
                                          check.node.labels=check.node.labels,
                                          ...)
                }
                else {
                    trees[[i]] <- phylo4(intreesphylolist[[i]],
                                         check.node.labels=check.node.labels,
                                         ...)
                }
            }
        }
        else {
            if (identical(check.node.labels, "asdata")) {
                if (is.null(intreesphylolist[[1]]$node.label)) {
                    warning("Could not use value \"asdata\" for ",
                            "check.node.labels because there are no ",
                            "labels associated with the tree ", i)
                    check.node.labels <- "drop"
                }
                trees <- phylo4d(intreesphylolist[[1]],
                                 check.node.labels=check.node.labels,
                                 ...)
            }
            else {
                trees <- phylo4(intreesphylolist[[1]],
                                check.node.labels=check.node.labels,
                                ...)
            }
        }
    }
    if (type == "tree" || (type == "all" && length(tipdata) == 0 )) {
        output <- trees
    }
    else {
        if (type == "data") {
            output <- tipdata
        }
        else {
            if (length(intreesphylolist) > 1 && !simplify) {
                output <- list()
                for (i in 1:length(intreesphylolist)) {
                    output[[i]] <- phylo4d(intreesphylolist[[i]],
                                           tip.data = tipdata,
                                           check.node.labels=check.node.labels,
                                           ...)
                }
            }
            else {
                output <- phylo4d(intreesphylolist[[1]],
                                  tip.data=tipdata,
                                  check.node.labels=check.node.labels,
                                  ...)
            }
        }
    }
    output
}

read.nexustreestring <- function(X) {
    ## Returns list of phylo objects (not multi.phylo, and always a list, even if
    ## there is only one element X is a character vector, each element is one line
    ## from a treefile
    ## This is based almost entirely on read.nexus from APE (Emmanuel Paradis).

    X<-unlist(strsplit(unlist(X),c("\n")))

    ## first remove all the comments

    ## BCO took out the "speedier removal of comments" code -- it keeps [&R]
    ## as a node label, replaced it with original APE code
    ## speedier removal of comments pc 13 April 2008
    ##X <- lapply(X, gsub, pattern = "\\[[^\\]]*\\]", replacement = "")

    LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)
    if (length(LEFT)) { # in case there are no comments at all
        w <- LEFT == RIGHT
        if (any(w)) { # in case all comments use at least 2 lines
            s <- LEFT[w]
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
            ## The above regexp was quite tough to find: it makes
            ## possible to delete series of comments on the same line:
            ##       ...[...]xxx[...]...
            ## without deleting the "xxx". This regexp is in three parts:
            ##       \\[      [^]]*       \\]
            ## where [^]]* means "any character, except "]", repeated zero
            ## or more times" (note that the ']' is not escaped here).
            ## The previous version was:
            ##       X[s] <- gsub("\\[.*\\]", "", X[s])
            ## which deleted the "xxx". (EP  2008-06-24)
        }
        w <- !w
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1))
            X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        }
    }
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- if (length(i2) == 1 && i2 > i1) TRUE else FALSE
    if (translation) {
        end <- semico[semico > i2][1]
        x <- X[(i2 + 1):end] # assumes there's a 'new line' after "TRANSLATE"
        ## x <- gsub("TRANSLATE", "", x, ignore.case = TRUE)
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
        n <- dim(TRANS)[1]
    }
    start <-
    if (translation)
        semico[semico > i2][1] + 1
    else
        semico[semico > i1][1]
    end <- endblock[endblock > i1][1] - 1
    tree <- X[start:end]
    rm(X)
    tree <- gsub("^.*= *", "", tree)
    semico <- grep(";", tree)
    Ntree <- length(semico)

    ## are some trees on several lines?
    if (any(diff(semico) != 1)) {
        STRING <- character(Ntree)
        s <- c(1, semico[-Ntree] + 1)
        j <- mapply(":", s, semico)
        for (i in 1:Ntree)
            STRING[i] <- paste(tree[j[, i]], collapse = "")
    }
    else
        STRING <- tree
    rm(tree)
    STRING <- gsub(" ", "", STRING)
    colon <- grep(":", STRING)
    if (!length(colon)) {
        ## TODO: recode clado.build, tree.build & .treeBuildWithTokens from ape to phylobase
        trees <- lapply(STRING, clado.build)
    }
    else {
        if (length(colon) == Ntree) {
            if (translation)
                trees <- lapply(STRING, .treeBuildWithTokens)
            else
               trees <- lapply(STRING, tree.build)

        }
        else {
            trees <- vector("list", Ntree)
            trees[colon] <- lapply(STRING[colon], tree.build)
            nocolon <- (1:Ntree)[!1:Ntree %in% colon]
            trees[nocolon] <- lapply(STRING[nocolon], clado.build)
            if (translation) {
                for (i in 1:Ntree) {
                    tr <- trees[[i]]
                    for (j in 1:n) {
                        ind <- which(tr$tip.label[j] == TRANS[, 1])
                        tr$tip.label[j] <- TRANS[ind, 2]
                    }
                    if (!is.null(tr$node.label)) {
                        for (j in 1:length(tr$node.label)) {
                            ind <- which(tr$node.label[j] == TRANS[, 1])
                            tr$node.label[j] <- TRANS[ind, 2]
                        }
                    }
                    trees[[i]] <- tr
                }
                translation <- FALSE
            }
        }
    }
    for (i in 1:Ntree) {
        tr <- trees[[i]]
        ## Check here that the root edge is not incorrectly represented
        ## in the object of class "phylo" by simply checking that there
        ## is a bifurcation at the root
        if (!translation) n <- length(tr$tip.label)
        ROOT <- n + 1
        if (sum(tr$edge[, 1] == ROOT) == 1 && dim(tr$edge)[1] > 1) {
            stop(paste("There is apparently two root edges in your file: ",
                       "cannot read tree file.\n  Reading NEXUS file aborted ",
                       "at tree no.", i, sep = ""))
            }
        }
    trees
}
