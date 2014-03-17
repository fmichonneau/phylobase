
R version 3.0.3 (2014-03-06) -- "Warm Puppy"
Copyright (C) 2014 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> ## RUnit script obtained from:
> ##   http://wiki.r-project.org/rwiki/doku.php?id=developers:runit
> 
> ## unit tests will not be done if RUnit is not available
> if(require("RUnit", quietly=TRUE)) {
+  
+   ## --- Setup ---
+  
+   pkg <- "phylobase"
+   if(Sys.getenv("RCMDCHECK") == "FALSE") {
+     ## Path to unit tests for standalone running under Makefile (not R CMD check)
+     ## PKG/tests/../inst/unitTests
+     path <- file.path(getwd(), "..", "inst", "unitTests")
+   } else {
+     ## Path to unit tests for R CMD check
+     ## PKG.Rcheck/tests/../PKG/unitTests
+     path <- system.file(package=pkg, "unitTests")
+   }
+   cat("\nRunning unit tests\n")
+   print(list(pkg=pkg, getwd=getwd(), pathToUnitTests=path))
+  
+   library(package=pkg, character.only=TRUE)
+  
+   ## If desired, load the name space to allow testing of private functions
+   ## if (is.element(pkg, loadedNamespaces()))
+   ##     attach(loadNamespace(pkg), name=paste("namespace", pkg, sep=":"), pos=3)
+   ##
+   ## or simply call PKG:::myPrivateFunction() in tests
+  
+   ## --- Testing ---
+  
+   ## Define tests
+   testSuite <- defineTestSuite(name=paste(pkg, "unit testing"),
+                                           dirs=path)
+   ## Run
+   tests <- runTestSuite(testSuite)
+  
+   ## Default report name
+   pathReport <- file.path(path, "report")
+  
+   ## Report to stdout and text files
+   cat("------------------- UNIT TEST SUMMARY ---------------------\n\n")
+   printTextProtocol(tests, showDetails=FALSE)
+   printTextProtocol(tests, showDetails=FALSE,
+                     fileName=paste(pathReport, "Summary.txt", sep=""))
+   printTextProtocol(tests, showDetails=TRUE,
+                     fileName=paste(pathReport, ".txt", sep=""))
+  
+   ## Report to HTML file
+   printHTMLProtocol(tests, fileName=paste(pathReport, ".html", sep=""))
+  
+   ## Return stop() to cause R CMD check stop in case of
+   ##  - failures i.e. FALSE to unit tests or
+   ##  - errors i.e. R errors
+   tmp <- getErrors(tests)
+   if(tmp$nFail > 0 | tmp$nErr > 0) {
+     stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
+                ", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
+   }
+ } else {
+   warning("cannot run unit tests -- package RUnit is not available")
+ }

Running unit tests
$pkg
[1] "phylobase"

$getwd
[1] "/home/francois/R-dev/phylobase/pkg/tests"

$pathToUnitTests
[1] "/home/francois/.R/library/phylobase/unitTests"


Warning messages:
1: replacing previous import by ‘Rcpp::evalCpp’ when loading ‘phylobase’ 
2: replacing previous import by ‘ade4::newick2phylog’ when loading ‘phylobase’ 
Execution halted
