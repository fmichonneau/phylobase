
[![Build Status](https://travis-ci.org/fmichonneau/phylobase.png?branch=master)](https://travis-ci.org/fmichonneau/phylobase) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/fmichonneau/phylobase?branch=master&svg=true)](https://ci.appveyor.com/project/fmichonneau/phylobase) [![codecov.io](http://codecov.io/github/fmichonneau/phylobase/coverage.svg?branch=master)](http://codecov.io/github/fmichonneau/phylobase?branch=master) ![](http://cranlogs.r-pkg.org/badges/phylobase) [![CRAN version](http://www.r-pkg.org/badges/version/phylobase)](http://cran.r-project.org/package=phylobase) [![Research software impact](http://depsy.org/api/package/cran/phylobase/badge.svg)](http://depsy.org/package/r/phylobase)

phylobase
=========

About this package
------------------

`phylobase` provides classes and methods to easily associate, manipulate, explore, and plot phylogenetic trees and data about the species they include. The goal of this package is to provide a base set of tools likely to be shared by all packages designed for phylogenetic analysis. This standardization will benefit both *end-users* by allowing them to move results across packages and keep data associated with the phylogenetic trees; and *developers* by focusing on method development instead of having to rewrite the base functions.

-   Authors: R Hackathon et al. (alphabetically: Ben Bolker, Marguerite Butler, Peter Cowan, Damien de Vienne, Dirk Eddelbuettel, Mark Holder, Thibaut Jombart, Steve Kembel, Francois Michonneau, David Orme, Brian O'Meara, Emmanuel Paradis, Jim Regetz, Derrick Zwickl)
-   Maintainer: Francois Michonneau
-   Licence: GPL (&gt;= 2)
-   Issues, bug reports, feature requests, discussion: <http://github.com/fmichonneau/phylobase/issues>

Installation
------------

### Stable version

The stable version (the minor and patch version numbers are even, e.g., 0.6.8) can be downloaded from CRAN.

``` r
install.packages("phylobase")
library(phylobase)
```

### Development version

The development version (the patch version number is odd, e.g., 0.6.9) is available on GitHub (<https://github.com/fmichonneau/phylobase>), and can be installed using the [`devtools`](http://cran.r-project.org/package=devtools) package.

``` r
library(devtools)
install_github("phylobase", "fmichonneau")
library(phylobase)
```

### Getting started

`phylobase` comes with example data sets `geospiza` and `geospiza_raw`.

-   `geospiza` is a `phylo4d` object (the `phylobase` class that holds together a phylogenetic tree and associated data, the `phylo4` class is for phylogenetic trees only).
-   `geospiza_raw` is a list that contains the tree `geospiza_raw$tree` (as an `ape::phylo` object) and the data `geospiza_raw$data` (as a `data.frame`) that were used to build the `geospiza` object.

Now we'll take the data from and merge it with the tree. However, since is included in the tree but not in the data set, we will initially run into some trouble:

``` r
data(geospiza_raw)
g1 <- as(geospiza_raw$tree, "phylo4")
geodata <- geospiza_raw$data
g2 <- phylo4d(g1, geodata)
```

    ## Error in formatData(phy = x, dt = tip.data, type = "tip", ...): The following nodes are not found in the dataset:  olivacea

To deal with *G. olivacea* missing from the data, we have a few choices. The easiest is to use `missing.data="warn"` to allow `R` to create the new object with a warning (you can also use `missing.data="OK"` to proceed without warnings):

``` r
g2 <- phylo4d(g1, geodata, missing.data="warn")
```

    ## Warning in formatData(phy = x, dt = tip.data, type = "tip", ...): The
    ## following nodes are not found in the dataset: olivacea

``` r
head(g2)
```

    ##           label node ancestor edge.length node.type    wingL  tarsusL
    ## 1    fuliginosa    1       24     0.05500       tip 4.132957 2.806514
    ## 2        fortis    2       24     0.05500       tip 4.244008 2.894717
    ## 3  magnirostris    3       23     0.11000       tip 4.404200 3.038950
    ## 4   conirostris    4       22     0.18333       tip 4.349867 2.984200
    ## 5      scandens    5       21     0.19250       tip 4.261222 2.929033
    ## 6    difficilis    6       20     0.22800       tip 4.224067 2.898917
    ## 7       pallida    7       25     0.08667       tip 4.265425 3.089450
    ## 8      parvulus    8       27     0.02000       tip 4.131600 2.973060
    ## 9    psittacula    9       27     0.02000       tip 4.235020 3.049120
    ## 10       pauper   10       26     0.03500       tip 4.232500 3.035900
    ##     culmenL    beakD   gonysW
    ## 1  2.094971 1.941157 1.845379
    ## 2  2.407025 2.362658 2.221867
    ## 3  2.724667 2.823767 2.675983
    ## 4  2.654400 2.513800 2.360167
    ## 5  2.621789 2.144700 2.036944
    ## 6  2.277183 2.011100 1.929983
    ## 7  2.430250 2.016350 1.949125
    ## 8  1.974420 1.873540 1.813340
    ## 9  2.259640 2.230040 2.073940
    ## 10 2.187000 2.073400 1.962100
    ##  [ reached getOption("max.print") -- omitted 10 rows ]

### Importing data

#### From NEXUS files

`phylobase` has a robust parser for NEXUS files (it uses the NEXUS Class Library from Paul Lewis and Mark Holder, [NCL](http://sourceforge.net/projects/ncl/files/)). It can be used to import simultaneously tree and species data.

``` r
myrmeFile <- system.file("nexusfiles/treeWithDiscreteData.nex", package="phylobase")
myrme <- readNexus(file=myrmeFile)
head(myrme)
```

    ##                      label node ancestor edge.length node.type        time
    ## 1   Myrmecocystussemirufus    1       27    1.724765       tip     diurnal
    ## 2   Myrmecocystusplacodops    2       27    1.724765       tip     diurnal
    ## 3      Myrmecocystusmendax    3       26    4.650818       tip     diurnal
    ## 4    Myrmecocystuskathjuli    4       28    1.083870       tip     diurnal
    ## 5    Myrmecocystuswheeleri    5       28    1.083870       tip     diurnal
    ## 6     Myrmecocystusmimicus    6       30    2.708942       tip     diurnal
    ## 7     Myrmecocystusdepilis    7       30    2.708942       tip     diurnal
    ## 8    Myrmecocystusromainei    8       32    2.193845       tip     diurnal
    ## 9  Myrmecocystusnequazcatl    9       32    2.193845       tip     diurnal
    ## 10       Myrmecocystusyuma   10       31    4.451425       tip crepuscular
    ## 11   Myrmecocystuskennedyi   11       23    6.044804       tip     diurnal
    ## 12 Myrmecocystuscreightoni   12       22   10.569191       tip crepuscular
    ## 13  Myrmecocystussnellingi   13       33    2.770378       tip crepuscular
    ## 14 Myrmecocystustenuinodis   14       33    2.770378       tip crepuscular
    ##         subgenus
    ## 1   Endiodioctes
    ## 2   Endiodioctes
    ## 3   Endiodioctes
    ## 4   Endiodioctes
    ## 5   Endiodioctes
    ## 6   Endiodioctes
    ## 7   Endiodioctes
    ## 8   Endiodioctes
    ## 9   Endiodioctes
    ## 10  Eremnocystus
    ## 11  Endiodioctes
    ## 12  Eremnocystus
    ## 13  Eremnocystus
    ## 14  Eremnocystus
    ##  [ reached getOption("max.print") -- omitted 6 rows ]

#### From NeXML

``` r
library(RNeXML)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:phylobase':
    ## 
    ##     edges

``` r
nxmlFile <- system.file("nexmlfiles/comp_analysis.xml", package="phylobase")
nxml <- nexml_read(nxmlFile)
nxmlEx <- phylo4(nxml)
```
