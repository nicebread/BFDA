<!-- README.md is generated from README.Rmd. Please edit that file 
rmarkdown::render("README.rmd")
-->

[![Build
Status](https://travis-ci.org/nicebread/BFDA.svg?branch=master)](https://travis-ci.org/nicebread/BFDA)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--02--17-yellowgreen.svg)](/commits/master)

[![packageversion](https://img.shields.io/badge/Package%20version-0.5.0-orange.svg?style=flat-square)](commits/master)

<!--
[![](https://codecov.io/gh/nicebread/BFDA/branch/master/graph/badge.svg)](https://codecov.io/gh/nicebread/BFDA)
-->

# BFDA Bayes factor design analysis

### Installation

The BFDA package is not on CRAN yet, but you can install the development
version from Github:

    library(devtools)
    install_github("nicebread/BFDA", subdir="package")

For installations on Windows the package requires R version 3.3.1 or
higher.

## How to use the BFDA package?

### 1\. Read our papers:

  - Schönbrodt, F. D. & Wagenmakers, E.-J. (2018). Bayes Factor Design
    Analysis: Planning for compelling evidence. *Psychonomic Bulletin &
    Review*, 25, 128-142. <doi:10.3758/s13423-017-1230-y>.
    \[[PDF](https://osf.io/d4dcu)\]\[[OSF project with reproducible
    code](https://osf.io/v7yxp/)\]
  - Stefan, A. M., Gronau, Q. F., Schönbrodt, F. D., & Wagenmakers, E.
    (2018). A Tutorial on Bayes Factor Design Analysis with Informed
    Priors. [PsyArXiv Preprint](https://doi.org/10.31234/osf.io/aqr79)

If you use this package to compute and report your design analysis,
please cite it as:

  - Schönbrodt, F. D. & Stefan, A. M. (2018). BFDA: An R package for
    Bayes factor design analysis (version 0.4.0). Retrieved from
    <https://github.com/nicebread/BFDA>

### 2\. Read the additional [manual](https://rawgit.com/nicebread/BFDA/master/package/doc/BFDA_manual.html).

![Sequential
Design](https://github.com/nicebread/BFDA/blob/master/movies/GIF1/BFDA1.gif)

## [BFDA in practice: A list of published examples](BFDA_examples.md)

## Extending the BFDA package

If you want to implement a new test (e.g., a BFDA for regression or
ANOVAs): The BFDA package uses a rather modular system for creating new
tests. For new tests, you need to:

1.  Add a new file with 4 functions:

<!-- end list -->

  - `sample.function`
  - `select.function`
  - `BF.test.function`
  - `freq.test.function`

Replace “function” with the test name, e.g.: sample.ANOVA etc.

2.  Add the new “type” (e.g., “ANOVA”) to the `print.BFDA` function in
    `R/1-Simulation.R`
3.  Add the new “type” (e.g., “ANOVA”) to the `BFDA.sanityCheck`
    function in `R/10-SanityChecks.R`
4.  Do a lot of testing\!

Probably it is easiest to take one of the existing test implementations
(e.g., `R/t.between.R` or `R/correlation.R`) and replace the name and
the content of the functions.

If you implement a new test, please let us know. The preferred workflow
would be that you fork the Github project, implement the test, do a lot
of testing, and send us a pull request. If everything works as intended,
we can add the test to the main project.
