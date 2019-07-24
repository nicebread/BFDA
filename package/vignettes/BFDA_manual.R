## ----setup, include=FALSE------------------------------------------------
	library(tint)
	knitr::opts_chunk$set(cache=TRUE, warnings = FALSE, messages=FALSE)
	# load all functions
	devtools::load_all()

## ----paired.t, warning=FALSE, fig.width=7, fig.height=4------------------
#devtools::install_github("nicebread/BFDA", subdir="package")
#library(BFDA)

# do a sequential design analysis
s1 <- BFDA.sim(expected.ES=0.4,
               prior=list("t", list(prior.location=0, prior.scale=sqrt(2)/2, prior.df=1)),
               n.min=50, stepsize=5, n.max=300, type="t.paired", design="sequential",
               alternative="greater", B=1000, cores=4, verbose=FALSE)
s0 <- BFDA.sim(expected.ES=0,
               prior=list("t", list(prior.location=0, prior.scale=sqrt(2)/2, prior.df=1)),
               n.min=50, stepsize=5, n.max=300, type="t.paired", design="sequential",
               alternative="greater", B=1000, cores=4, verbose=FALSE)

# if no n.min and n.max is provided in the `BFDA.analyze` function,
# the values from the simulation are taken
BFDA.analyze(s1, design="sequential", boundary=10)
BFDA.analyze(s0, design="sequential", boundary=10)

BFDA.analyze(s1, design="sequential", boundary=6)
BFDA.analyze(s0, design="sequential", boundary=6)

plot(s1)

## ----corr_walkthrough, warning=FALSE, fig.width=7, fig.height=4----------
#devtools::install_github("nicebread/BFDA", subdir="package")
#library(BFDA)

# do a sequential design analysis
c1 <- BFDA.sim(expected.ES=0.21, prior=list("stretchedbeta",list(prior.kappa=2)),
               n.min=50, stepsize=10, n.max=300, B=1000, type="correlation",
               design="sequential", alternative="greater", cores=4, verbose=FALSE)
c0 <- BFDA.sim(expected.ES=0, prior=list("stretchedbeta", list(prior.kappa=2)),
               n.min=50, stepsize=10, n.max=300, B=1000, type="correlation",
               design="sequential", alternative="greater", cores=4, verbose=FALSE)

# if no n.min and n.max is provided in the `BFDA.analyze` function,
# the values from the simulation are taken
BFDA.analyze(c1, design="sequential", boundary=10)
BFDA.analyze(c0, design="sequential", boundary=10)

BFDA.analyze(c1, design="sequential", boundary=6)
BFDA.analyze(c0, design="sequential", boundary=6)

plot(c1, boundary=c(1/10, 20), n.max=150)

## ----abtest_walkthrough, warning=FALSE, fig.width=7, fig.height=4--------
#devtools::install_github("nicebread/BFDA", subdir="package")
#library(BFDA)

# do a sequential design analysis
ab1 <- BFDA.sim(expected.ES=2, prior=list("normal", list(prior.mean=0.5, prior.variance = 1)),
                n.min=50, stepsize=10, n.max=300, B=1000, type="abtest",
                design="sequential", alternative="two.sided", cores=4, verbose=FALSE,
                options.sample = list(effecttype = "OR"))
ab0 <- BFDA.sim(expected.ES=0, prior=list("normal", list(prior.mean=0.5, prior.variance = 1)),
                n.min=50, stepsize=10, n.max=300, B=1000, type="abtest",
                design="sequential", alternative="two.sided", cores=4, verbose=FALSE,
                options.sample = list(effecttype = "OR"))

# if no n.min and n.max is provided in the `BFDA.analyze` function,
# the values from the simulation are taken
BFDA.analyze(ab1, design="sequential", boundary=10)
BFDA.analyze(ab0, design="sequential", boundary=10)

BFDA.analyze(ab1, design="sequential", boundary=6)
BFDA.analyze(ab0, design="sequential", boundary=6)

plot(ab1, boundary=c(1/10, 20), n.max=150)

