[![Build Status](https://travis-ci.org/nicebread/BFDA.svg?branch=master)](https://travis-ci.org/nicebread/BFDA)

# BFDA Bayes factor design analysis #

### Installation

The BFDA package is not on CRAN yet, but you can install the development version from Github:

    library(devtools)
    install_github("nicebread/BFDA", subdir="package")

For installations on Windows the package requires R version 3.3.1 or higher.

### How to use the BFDA package?

1. Read our [paper](http://papers.ssrn.com/abstract=2722435):

	> Schönbrodt, F. D. & Wagenmakers, E.-J. (submitted). Bayes Factor Design Analysis: Planning for compelling evidence. Retrieved from http://ssrn.com/abstract=2722435.

2. Read the additional [manual](https://rawgit.com/nicebread/BFDA/master/vignette/BFDA_manual.html).

![Sequential Design](https://github.com/nicebread/BFDA/blob/master/movies/GIF1/BFDA1.gif)


## Extending the BFDA package

If you want to implement a new test (e.g., a BFDA for regression or ANOVAs): The BFDA package uses a rather modular system for creating new tests. For new tests, you need to:

1. Add a new file with 4 functions:

- `sample.function`
- `select.function`
- `BF.test.function`
- `freq.test.function`

Replace "function" with the test name, e.g.: sample.ANOVA etc.

2. Add the new "type" (e.g., "ANOVA") to the `print.BFDA` function in `R/1-Simulation.R`
3. Add the new "type" (e.g., "ANOVA") to the `BFDA.sanityCheck` function in `R/10-SanityChecks.R`
4. Do a lot of testing!

Probably it is easiest to take one of the existing test implementations (e.g., `R/t.between.R` or `R/correlation.R`)  and replace the name and the content of the functions.

If you implement a new test, please let me know. The preferred workflow would be that you fork the Github project, implement the test, and send me a pull request. If everything works as intended, we can add the test to the main project.