phylobase.options <- function (...) {
    ## code from package sm
    if (nargs() == 0) return(.phylobase.Options)
    current <- .phylobase.Options
    if (is.character(...))
        temp <- eval(parse(text = paste(c("list(", ..., ")"))))
    else temp <- list(...)
    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
               list = temp <- arg,
               character = return(.phylobase.Options[arg]),
               stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")
    changed <- current[n]
    current[n] <- temp
    if (sys.parent() == 0) env <- asNamespace("phylobase") else env <- parent.frame()
    assign(".phylobase.Options", current, envir = env)
    invisible(current)
}
