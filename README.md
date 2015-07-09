

[![Build Status](https://travis-ci.org/fmichonneau/phylobase.png?branch=master)](https://travis-ci.org/fmichonneau/phylobase)
[![Coverage Status](https://coveralls.io/repos/fmichonneau/phylobase/badge.svg?branch=master&service=github)](https://coveralls.io/github/fmichonneau/phylobase?branch=master)
![](http://cranlogs.r-pkg.org/badges/phylobase)
[![CRAN version](http://www.r-pkg.org/badges/version/phylobase)](http://cran.r-project.org/package=phylobase)

# phylobase

## About this package

`phylobase` provides classes and methods to easily associate, manipulate,
explore, and plot phylogenetic trees and data about the species they
include. The goal of this package is to provide a base set of tools likely to be
shared by all packages designed for phylogenetic analysis. This standardization
will benefit both *end-users* by allowing them to move results across packages
and keep data associated with the phylogenetic trees; and *developers* by
focusing on method development instead of having to rewrite the base functions.

- Authors: R Hackathon et al. (alphabetically: Ben Bolker, Marguerite Butler,
  Peter Cowan, Damien de Vienne, Dirk Eddelbuettel, Mark Holder, Thibaut
  Jombart, Steve Kembel, Francois Michonneau, David Orme, Brian O'Meara,
  Emmanuel Paradis, Jim Regetz, Derrick Zwickl)
- Maintainer: Francois Michonneau
- Licence: GPL (>= 2)
- Issues, bug reports, feature requests, discussion:
  http://github.com/fmichonneau/phylobase/issues

## Installation

### Stable version

The stable version (the minor and patch version numbers are even, e.g., 0.6.8)
can be downloaded from CRAN.


```r
install.packages("phylobase")
library(phylobase)
```

### Development version

The development version (the patch version number is odd, e.g., 0.6.9) is
available on GitHub (https://github.com/fmichonneau/phylobase), and can be
installed using the [`devtools`](http://cran.r-project.org/package=devtools)
package.


```r
library(devtools)
install_github("phylobase", "fmichonneau")
library(phylobase)
```

### Getting started



`phylobase` comes with example data sets `geospiza` and `geospiza_raw`.

- `geospiza` is a `phylo4d` object (the `phylobase` class that holds together a
  phylogenetic tree and associated data, the `phylo4` class is for phylogenetic
  trees only).
- `geospiza_raw` is a list that contains the tree `geospiza_raw$tree` (as an
  `ape::phylo` object) and the data `geospiza_raw$data` (as a `data.frame`) that
  were used to build the `geospiza` object.

Now we'll take the \emph{Geospiza} data from \verb+geospiza_raw$data+ and merge
it with the tree. However, since \emph{G. olivacea} is included in the tree but
not in the data set, we will initially run into some trouble:


```r
data(geospiza_raw)
g1 <- as(geospiza_raw$tree, "phylo4")
geodata <- geospiza_raw$data
g2 <- phylo4d(g1, geodata)
```

```
## Error in formatData(phy = x, dt = tip.data, type = "tip", ...): The following nodes are not found in the dataset:  olivacea
```

To deal with _G. olivacea_ missing from the data, we have a few choices. The
easiest is to use `missing.data="warn"` to allow `R` to create the new object
with a warning (you can also use `missing.data="OK"` to proceed without
warnings):


```r
g2 <- phylo4d(g1, geodata, missing.data="warn")
```

```
## Warning in formatData(phy = x, dt = tip.data, type = "tip", ...): The
## following nodes are not found in the dataset: olivacea
```

```r
head(g2)
```

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
## 11   Platyspiza   11       18     0.46550       tip 4.419686 3.270543
## 12        fusca   12       17     0.53409       tip 3.975393 2.936536
## 13 Pinaroloxias   13       16     0.58333       tip 4.188600 2.980200
## 14     olivacea   14       15     0.88077       tip       NA       NA
## 15         <NA>   15        0          NA      root       NA       NA
## 16         <NA>   16       15     0.29744  internal       NA       NA
## 17         <NA>   17       16     0.04924  internal       NA       NA
## 18         <NA>   18       17     0.06859  internal       NA       NA
## 19         <NA>   19       18     0.13404  internal       NA       NA
## 20         <NA>   20       19     0.10346  internal       NA       NA
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
## 11 2.331471 2.347471 2.282443
## 12 2.051843 1.191264 1.401186
## 13 2.311100 1.547500 1.630100
## 14       NA       NA       NA
## 15       NA       NA       NA
## 16       NA       NA       NA
## 17       NA       NA       NA
## 18       NA       NA       NA
## 19       NA       NA       NA
## 20       NA       NA       NA
```

### Importing data

#### From NEXUS files

`phylobase` has a robust parser for NEXUS files (it uses the NEXUS Class Library
from Paul Lewis and Mark Holder,
[NCL](http://sourceforge.net/projects/ncl/files/)). It can be used to import
simultaneously tree and species data.


```r
myrmeFile <- system.file("nexusfiles/treeWithDiscreteData.nex", package="phylobase")
myrme <- readNexus(file=myrmeFile)
head(myrme)
```

```
##                      label node ancestor edge.length node.type        time
## 1    Myrmecocystuscfnavajo    1       35    2.869547       tip   nocturnal
## 2  Myrmecocystuscreightoni    2       22   10.569191       tip crepuscular
## 3     Myrmecocystusdepilis    3       30    2.708942       tip     diurnal
## 4    Myrmecocystuskathjuli    4       28    1.083870       tip     diurnal
## 5    Myrmecocystuskennedyi    5       23    6.044804       tip     diurnal
## 6      Myrmecocystusmendax    6       26    4.650818       tip     diurnal
## 7   Myrmecocystusmexicanus    7       34    5.724923       tip   nocturnal
## 8     Myrmecocystusmimicus    8       30    2.708942       tip     diurnal
## 9      Myrmecocystusnavajo    9       35    2.869547       tip   nocturnal
## 10 Myrmecocystusnequazcatl   10       32    2.193845       tip     diurnal
## 11  Myrmecocystusplacodops   11       27    1.724765       tip     diurnal
## 12   Myrmecocystusromainei   12       32    2.193845       tip     diurnal
## 13  Myrmecocystussemirufus   13       27    1.724765       tip     diurnal
## 14  Myrmecocystussnellingi   14       33    2.770378       tip crepuscular
## 15 Myrmecocystustenuinodis   15       33    2.770378       tip crepuscular
## 16  Myrmecocystustestaceus   16       20   12.300701       tip crepuscular
## 17   Myrmecocystuswheeleri   17       28    1.083870       tip     diurnal
## 18       Myrmecocystusyuma   18       31    4.451425       tip crepuscular
## 19                    <NA>   19        0          NA      root        <NA>
## 20                    <NA>   20       19    1.699299  internal        <NA>
##         subgenus
## 1  Myrmecocystus
## 2   Eremnocystus
## 3   Endiodioctes
## 4   Endiodioctes
## 5   Endiodioctes
## 6   Endiodioctes
## 7  Myrmecocystus
## 8   Endiodioctes
## 9  Myrmecocystus
## 10  Endiodioctes
## 11  Endiodioctes
## 12  Endiodioctes
## 13  Endiodioctes
## 14  Eremnocystus
## 15  Eremnocystus
## 16 Myrmecocystus
## 17  Endiodioctes
## 18  Eremnocystus
## 19          <NA>
## 20          <NA>
```

#### From NeXML


```r
library(RNeXML)
```

```
## Loading required package: ape
## 
## Attaching package: 'ape'
## 
## The following object is masked from 'package:phylobase':
## 
##     edges
```

```r
nxmlFile <- system.file("nexmlfiles/comp_analysis.xml", package="phylobase")
nxml <- nexml_read(nxmlFile)
nxmlEx <- phylo4(nxml)
```
