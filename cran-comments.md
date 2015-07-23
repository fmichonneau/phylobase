## Test environments

- local Ubuntu 15.04, R 3.2.1
- Ubuntu 12.04 (travis-ci), R 3.2.1
- Windows with win-builder (R 3.2.1 and R-devel r68715)
- Debian Testing with R-devel (r68728) compiled with gcc-5

## R CMD check results

- There were no ERRORs or WARNINGs

- There was 1 NOTE:
  * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Francois Michonneau <francois.michonneau@gmail.com>'
    Possibly mis-spelled words in DESCRIPTION:
    Phylogenetic (4:25)

## Downstream dependencies

I have run R CMD check on all downstream dependencies listed on CRAN. All
packages passed. The only package that produced a WARNING relevant to the
changes in phylobase is `phyloTop` as this package uses a function that is
deprecated in the present version. I contacted the maintainer of this package.
