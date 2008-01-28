ReadWithNCLVersion <- function() {

  version <- installed.packages()["ReadWithNCL","Version"]
                                 
  # Print version of this package
  cat('\nReadWithNCL package Version ', version, '\n\n')

  # Print version of Rcpp/RcppTemplate used to build this package
  licenseFile <- file(system.file(".","LICENSE-Rcpp.txt",package="ReadWithNCL"),"r")
  writeLines(readLines(licenseFile)[1:4])
}
