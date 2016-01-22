This submission follows an email sent to me by Prof. Brian Ripley and Martin
Maelcher regarding a typo in the examples (detected by changes to R-devel).

## Test environments

- local Ubuntu 15.04, R 3.2.3
- Ubuntu 12.04 (travis-ci), R 3.2.3
- Windows with win-builder (R 3.2.3 and R-devel r69971)
- Debian Testing with R-devel compiled with gcc-5

## R CMD check results

- There were no ERRORs or WARNINGs

- There was 1 NOTE:
  * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Francois Michonneau <francois.michonneau@gmail.com>'
    Possibly mis-spelled words in DESCRIPTION:
    Phylogenetic (4:25)

## Downstream dependencies

I have run R CMD check on all downstream dependencies listed on CRAN. None of
the packages seem to have generated NOTEs or WARNINGs relevant to the
changes in phylobase.
