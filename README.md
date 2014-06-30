

[![Build Status](https://travis-ci.org/fmichonneau/phylobase.png?branch=master)](https://travis-ci.org/fmichonneau/phylobase.png)

Latest Build log: https://travis-ci.org/fmichonneau/phylobase

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
installed using the `[devtools](http://cran.r-project.org/package=devtools)`
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
## Error: The following nodes are not found in the dataset: olivacea
```

To deal with _G. olivacea_ missing from the data, we have a few choices. The
easiest is to use `missing.data="warn"` to allow `R` to create the new object
with a warning (you can also use `missing.data="OK"` to proceed without
warnings):


```r
g2 <- phylo4d(g1, geodata, missing.data="warn")
```

```
## Warning: The following nodes are not found in the dataset: olivacea
```

```r
head(g2)
```

```
##           label node ancestor edge.length node.type wingL tarsusL culmenL
## 1    fuliginosa    1       24     0.05500       tip 4.133   2.807   2.095
## 2        fortis    2       24     0.05500       tip 4.244   2.895   2.407
## 3  magnirostris    3       23     0.11000       tip 4.404   3.039   2.725
## 4   conirostris    4       22     0.18333       tip 4.350   2.984   2.654
## 5      scandens    5       21     0.19250       tip 4.261   2.929   2.622
## 6    difficilis    6       20     0.22800       tip 4.224   2.899   2.277
## 7       pallida    7       25     0.08667       tip 4.265   3.089   2.430
## 8      parvulus    8       27     0.02000       tip 4.132   2.973   1.974
## 9    psittacula    9       27     0.02000       tip 4.235   3.049   2.260
## 10       pauper   10       26     0.03500       tip 4.232   3.036   2.187
## 11   Platyspiza   11       18     0.46550       tip 4.420   3.271   2.331
## 12        fusca   12       17     0.53409       tip 3.975   2.937   2.052
## 13 Pinaroloxias   13       16     0.58333       tip 4.189   2.980   2.311
## 14     olivacea   14       15     0.88077       tip    NA      NA      NA
## 15         <NA>   15        0          NA      root    NA      NA      NA
## 16         <NA>   16       15     0.29744  internal    NA      NA      NA
## 17         <NA>   17       16     0.04924  internal    NA      NA      NA
## 18         <NA>   18       17     0.06859  internal    NA      NA      NA
## 19         <NA>   19       18     0.13404  internal    NA      NA      NA
## 20         <NA>   20       19     0.10346  internal    NA      NA      NA
##    beakD gonysW
## 1  1.941  1.845
## 2  2.363  2.222
## 3  2.824  2.676
## 4  2.514  2.360
## 5  2.145  2.037
## 6  2.011  1.930
## 7  2.016  1.949
## 8  1.874  1.813
## 9  2.230  2.074
## 10 2.073  1.962
## 11 2.347  2.282
## 12 1.191  1.401
## 13 1.548  1.630
## 14    NA     NA
## 15    NA     NA
## 16    NA     NA
## 17    NA     NA
## 18    NA     NA
## 19    NA     NA
## 20    NA     NA
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
## 1   Myrmecocystussemirufus    1       27       1.725       tip     diurnal
## 2   Myrmecocystusplacodops    2       27       1.725       tip     diurnal
## 3      Myrmecocystusmendax    3       26       4.651       tip     diurnal
## 4    Myrmecocystuskathjuli    4       28       1.084       tip     diurnal
## 5    Myrmecocystuswheeleri    5       28       1.084       tip     diurnal
## 6     Myrmecocystusmimicus    6       30       2.709       tip     diurnal
## 7     Myrmecocystusdepilis    7       30       2.709       tip     diurnal
## 8    Myrmecocystusromainei    8       32       2.194       tip     diurnal
## 9  Myrmecocystusnequazcatl    9       32       2.194       tip     diurnal
## 10       Myrmecocystusyuma   10       31       4.451       tip crepuscular
## 11   Myrmecocystuskennedyi   11       23       6.045       tip     diurnal
## 12 Myrmecocystuscreightoni   12       22      10.569       tip crepuscular
## 13  Myrmecocystussnellingi   13       33       2.770       tip crepuscular
## 14 Myrmecocystustenuinodis   14       33       2.770       tip crepuscular
## 15  Myrmecocystustestaceus   15       20      12.301       tip crepuscular
## 16  Myrmecocystusmexicanus   16       34       5.725       tip   nocturnal
## 17   Myrmecocystuscfnavajo   17       35       2.870       tip   nocturnal
## 18     Myrmecocystusnavajo   18       35       2.870       tip   nocturnal
## 19                    <NA>   19        0          NA      root        <NA>
## 20                    <NA>   20       19       1.699  internal        <NA>
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
## 15 Myrmecocystus
## 16 Myrmecocystus
## 17 Myrmecocystus
## 18 Myrmecocystus
## 19          <NA>
## 20          <NA>
```

#### From NeXML


```r
library(RNeXML)
nxmlFile <- system.file("nexmlfiles/comp_analysis.xml", package="phylobase")
nxml <- nexml_read(nxmlFile)
nxmlEx <- phylo4(nxml)
```
