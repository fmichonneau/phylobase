readNexus <- function (file, simplify=FALSE, type=c("all", "tree", "data"),
                       char.all=FALSE, polymorphic.convert=TRUE,
                       levels.uniform=TRUE, quiet=TRUE,
                       check.node.labels=c("keep", "drop", "asdata"),
                       return.labels=TRUE, ...) {

  return(readNCL(file=file, simplify=simplify, type=type, char.all=char.all,
          polymorphic.convert=polymorphic.convert, levels.uniform=levels.uniform,
          quiet=quiet, check.node.labels=check.node.labels,
          return.labels=return.labels, ...))
}
