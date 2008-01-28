# Init file for package RcppTemplate
.First.lib <- function(lib, pkg) {
	require(phylo4)
	require(ape)
  library.dynam("RcppTemplate", pkg, lib )
}

# Runs RcppTemplate demos interactively.
RcppTemplateDemo <- function() {
  demoCode <- demo(package="RcppTemplate")$results[,3]
  if(length(demoCode) == 0) {
    return('No demos found in package RcppTemplate')
  }
  while(TRUE) {
    for(i in seq(1,length(demoCode))) {
      cat(i, ': ', demoCode[i], '\n')
    }
    suppressWarnings(demoNum <- as.integer(readline(prompt="Enter demo number to run [0 to view descriptions, RETURN to exit]: ")))
    
    if(is.na(demoNum)) {
      cat("Didn't get a number, exiting.\n")
      break
    }
    else if(demoNum == 0) {
      print(demo(package="RcppTemplate"))
    }
    else {
      if(demoNum < 0 || demoNum > length(demoCode)) {
        cat('Number out of range, try again...\n')
        next
      }
      demo(demoCode[demoNum], character.only=TRUE)
    }
  }
}

