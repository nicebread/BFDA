BPA.d.0.5one <- BPA.sim.ttest(d=0.5, n.min=20, n.max=300, boundary=Inf, design="sequential", B=100, verbose=TRUE, cores=2, rscale=1, alternative="directional", ETA=FALSE)

BPA.d.0.5two <- BPA.sim.ttest(d=0.5, n.min=20, n.max=300, boundary=Inf, design="sequential", B=100, verbose=TRUE, cores=2, rscale=1, alternative="un", ETA=FALSE)


# ---------------------------------------------------------------------
# Simulate data
BPA.d.0.5 <- BPA.sim.ttest(d=0.5, n.min=20, n.max=500, boundary=Inf, stepsize=1, design="sequential", B=10000, verbose=TRUE, cores=2, rscale=1)
save(BPA.d.0.5, file="../finalSims/sim.0.5.RData")

BPA.d.0.5.prior <- BPA.sim.ttest(d=rnorm(100000, 0.5, 0.1), n.min=20, n.max=500, boundary=Inf, stepsize=1, design="sequential", B=10000, verbose=TRUE, cores=2, rscale=1)
save(BPA.d.0.5.prior, file="../finalSims/sim.0.5.prior.RData")

BPA.d.0 <- BPA.sim.ttest(d=0, n.min=20, n.max=500, boundary=Inf, stepsize=1, design="sequential", B=10000, verbose=TRUE, cores=2, rscale=1)
save(BPA.d.0, file="../finalSims/sim.0.RData")






plot(BPA.d.0.5, boundar=Inf, n.max=70)

load("../finalSims/sim0.5_step1.RData")
plot(sim, n.max=40, boundary=70, dens.amplification=1, xextension=1.5)

plot(sim, n.max=40, boundary=Inf)

ind <- BPA.analysis(sim, n.max=4)

plot(sim, n.max=35, boundary=4, dens.amplification=1)

plot(sim, n.max=60, boundary=4, dens.amplification=1.5, n.trajectories=0)


BPA.analysis(sim, n.max=60, boundary=4)
plotSim(sim, n.max=30, boundary=Inf, traj.selection="fixed", xlim=c(20, 400))

plotSim(sim, n.max=30, boundary=Inf, traj.selection="fixed", ylim=c(log(1/30), log(1010)))


## ======================================================================
## PLot the three plots for the paper
## ======================================================================

# Figure 1: fixed-n
load("../finalSims/sim0.5_step1.RData")
pdf(file="Figure1.pdf", width=7, height=4, pointsize=10)

plot(sim, n.max=60, boundary=Inf, traj.selection="prop", ylim=c(log(1/35), log(10000)), dens.amplification=0.5, xextension=1.8, yaxis.at=c(-log(30), -log(10), -log(3), log(1), log(3), log(10), log(30), log(100), log(1000), log(10000)), yaxis.labels=c("1/30", "1/10", "1/3", "1", "3", "10", "30", "100", "1000", "10000"))

dev.off()



# Figure 2: unlimited SBF
pdf(file="Figure2.pdf", width=7, height=4, pointsize=10)
plot(sim, n.max=NA, xlim=c(10, 240), boundary=6, traj.selection="fixed", xextension=1.8)
dev.off()


# Figure 3: SBF+maxN
pdf(file="Figure3.pdf", width=7, height=4, pointsize=10)
plot(sim, n.max=60, xlim=c(10, 240), boundary=10, traj.selection="fixed", xextension=0.8)
dev.off()



# ---------------------------------------------------------------------
# SSD plot test

#load("../../finalSims/sim0.5_step1.RData")

load("../finalSims/sim.0.5.prior.RData") #sim3

# Correct sign of BF?
SSD(sim3, power=.80, criterion=c(1/3, 3))


# symmetric boundaries at 1/3 and 3
SSD(sim, power=.90, criterion=c(1/3, 3))

# asymmetric boundaries at 0.5 and 10
SSD(sim, power=.80, criterion=c(0.5, 10))

# decrease Type-II errors with stronger
SSD(sim, power=.80, criterion=c(1/4, 10))


load("../finalSims/sim.0.RData")
load("../finalSims/sim.0.3.RData")
load("../finalSims/sim.0.5.prior.RData")

SSD(sim, power=.50, criterion=c(1/3, 3))

# Prepare lean ES=0 object
load("../finalSims/sim.0.RData")
BPA0 <- sim2
BPA0$sim <- BPA0$sim[, c("d", "boundary", "n", "logBF", "d.emp", "id")]
