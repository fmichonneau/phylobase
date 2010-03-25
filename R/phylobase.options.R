### Modified code from package sm
phylobase.options <- function (...) {
    if (nargs() == 0) return(.phylobase.Options)
    current <- .phylobase.Options
    temp <- list(...)
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

    if (!all(names(temp) %in% names(current)))
        stop("Option name invalid: ", sQuote(names(temp)))
    changed <- current[n]
    current[n] <- temp
    current <- lapply(current, function(foo) {
        foo <- match.arg(foo, c("warn", "fail", "ok"))
    })
    ## options are always global
    env <- asNamespace("phylobase")
    assign(".phylobase.Options", current, envir = env)
    invisible(current)
}
