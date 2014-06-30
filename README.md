

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
