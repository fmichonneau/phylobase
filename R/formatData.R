formatData <- function(phy, dt, type=c("tip", "internal", "all"),
                       match.data=TRUE, rownamesAsLabels=FALSE,
                       label.type=c("rownames", "column"),
                       label.column=1, missing.data=c("fail", "warn", "OK"),
                       extra.data=c("warn", "OK", "fail")
                       ) {

    type <- match.arg(type)
    label.type <- match.arg(label.type)
    stopifnot(label.column %in% 1:ncol(dt))
    missing.data <- match.arg(missing.data)
    extra.data <- match.arg(extra.data)

    nr <- switch(type,
                 tip = nTips(phy),
                 internal = nNodes(phy),
                 all = nTips(phy)+nNodes(phy))

    tmpDt <- array(, dim=c(nr, ncol(dt)),
                   dimnames=list(nodeId(phy, type), colnames(dt)))
    tmpDt <- data.frame(tmpDt)

    if(match.data) {
        ## Replace node labels by node numbers
        ndNames <- switch(label.type,
                          rownames = rownames(dt),
                          column = dt[,label.column])
        ndDt <- lapply(ndNames, function(nd) {
            if(nchar(gsub("[0-9]", "", nd)) == 0 && !rownamesAsLabels)
                getNode(phy, as.integer(nd), missing="OK")
            else getNode(phy, nd, missing="OK")
        })
        ndDt <- unlist(ndDt)

        ## Make sure that data are matched to appropriate nodes
        if(type != "all") {
            switch(type,
                   tip = {
                     ## BMB: don't bother trying to match NAs
                       if(any(na.omit(names(ndDt)) %in% labels(phy, "internal")))
                           stop("You are trying to match tip data to internal ",
                                "nodes. Make sure that your data identifiers ",
                                "are correct.")
                   },
                   internal = {
                       if(any(na.omit(names(ndDt)) %in% labels(phy, "tip")))
                           stop("You are trying to match node data to tip ",
                                "nodes. Make sure that your data identifiers ",
                                "are correct.")
                   })
        }

        ## Check differences
        extra <- names(ndDt[is.na(ndDt)])
        mssng <- nodeId(phy, type)[! nodeId(phy, type) %in% ndDt]

        if(length(mssng) > 0 && missing.data != "OK") {
            msg <- "The following nodes are not found in the dataset: "

            ## provides label if it exists and node number otherwise
            mssng <- sapply(mssng, function(m) {
                m <- getNode(phy, m)
                if (is.na(names(m)) || is.null(names(m)))
                    m
                else
                    names(m)
            })

            msg <- paste(msg, paste(mssng, collapse=", "))
            switch(missing.data,
                   warn = warning(msg),
                   fail = stop(msg))
        }

        if(length(extra) > 0 && extra.data != "OK") {
            msg <- "The following names are not found in the tree: "
            msg <- paste(msg, paste(extra, collapse=", "))
            switch(extra.data,
                   warn = warning(msg),
                   fail = stop(msg))

        }
        ## Format data to have correct dimensions
        dt <- dt[!is.na(ndDt) ,, drop=FALSE]
        rownames(dt) <- ndDt[!is.na(ndDt)]
        if(label.type == "column") dt <- dt[, -label.column]
        tmpDt[,] <- dt[match(rownames(tmpDt), rownames(dt)) ,, drop=FALSE]
    }
    else {
        ## Remove rownames in data provided
        rownames(dt) <- NULL

        ## Tips before internal nodes for all.data
        if (type == "all")
            rownames(tmpDt) <- 1:nr

        ## Check differences between dataset and tree
        diffNr <- nrow(dt) - nr
        if(diffNr > 0 && extra.data != "OK") {
            msg <- paste("There are", diffNr, "extra rows.")
            switch(extra.data,
                   warn = warning(msg),
                   fail = stop(msg))
        }
        if(diffNr < 0 && missing.data != "OK") {
            msg <- paste("There are", abs(diffNr), "missing rows.")
            switch(missing.data,
                   warn = warning(msg),
                   fail = stop(msg))
        }
        tmpDt[,] <- dt[1:min(nrow(dt), nr) ,, drop = FALSE]
    }

    tmpDt
}
