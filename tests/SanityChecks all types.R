

# sim <- BFDA.sim(d=0.5, n.min=20, n.max=300, boundary=Inf, stepsize=1, design="sequential", B=1000, verbose=TRUE, cores=2)
# save(sim, file="sim.0.5.RData")
#
# sim <- BFDA.ttest.2(d=rnorm(100000, 0.5, 0.1), n.min=20, n.max=1000, boundary=Inf, stepsize=NA, design="sequential", B=10000, verbose=TRUE, cores=2)
# save(sim, file="sim.prior.0.5.RData")
#
# sim <- BFDA.ttest.2(d=0, n.min=20, n.max=1000, boundary=Inf, stepsize=NA, design="sequential", B=10000, verbose=TRUE, cores=2)
# save(sim, file="sim.0.RData")

# sim <- BFDA.ttest.2(d=0.5, n.min=20, n.max=400, boundary=Inf, stepsize=2, design="sequential", B=5000, verbose=TRUE, cores=2)
# save(sim, file="sim.0.5.step2RData")


# ---------------------------------------------------------------------
#

#sim <- BFDA.sim(d=0, n.min=20, n.max=500, boundary=Inf, stepsize=NA, design="sequential", B=1000, verbose=TRUE, cores=2)
#save(sim, file="../../finalSims/sim.0b.RData")


# sim1 <- BFDA.sim(expected.ES=0.5, type="t.between", n.min=20, n.max=100, boundary=Inf, stepsize=NA, design="sequential", B=10, verbose=TRUE, cores=1, ETA=FALSE)

#sim0 <- BFDA.sim(expected.ES=0, type="between", n.min=20, n.max=100, boundary=Inf, stepsize=NA, design="sequential", B=10, verbose=TRUE, cores=1, ETA=FALSE)

#evDens(sim1, sim0, n=50)
#SSD(sim1, power=.80)

## ======================================================================
## Correlation
## ======================================================================


corr.H1 <- BFDA.sim(expected.ES=0.25, type="correlation", n.min=20, n.max=400, boundary=Inf, stepsize=10, design="sequential", B=1000, verbose=TRUE, cores=1, ETA=FALSE)
save(corr.H1, file="simfiles/corr.H1.RData")
load(file="simfiles/corr.H1.RData")

corr.H0 <- BFDA.sim(expected.ES=0, type="correlation", n.min=20, n.max=1000, boundary=Inf, stepsize=10, design="sequential", B=1000, verbose=TRUE, cores=1, ETA=FALSE)
save(corr.H0, file="simfiles/corr.H0.RData")
load(file="simfiles/corr.H0.RData")

plot(corr.H1, boundary=6)
BFDA.sanityCheck(corr.H1)
BFDA.analyze(corr.H1, n.min=50, n.max=50, boundary=3)

plot(corr.H0, boundary=6)
BFDA.sanityCheck(corr.H0)
BFDA.analyze(corr.H0, n.min=50, n.max=50, boundary=6)



# ---------------------------------------------------------------------
# check timings

corr.H0a <- BFDA.sim(expected.ES=0, type="correlation", n.min=20, n.max=1000, boundary=Inf, stepsize=10, design="sequential", B=1000, verbose=TRUE, cores=1, ETA=FALSE)

corr.H0b <- BFDA.sim(expected.ES=0, type="correlation", n.min=20, n.max=1000, boundary=Inf, stepsize=10, design="sequential", B=1000, verbose=TRUE, cores=4, ETA=FALSE)