


load("../finalSims/sim0.5_step1.RData")
plot(sim, n.max=60, boundary=4, dens.amplification=1.5)

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

plot(sim, n.max=60, boundary=Inf, traj.selection="prop", ylim=c(log(1/35), log(10000)), xextension=1.8, yaxis.at=c(-log(30), -log(10), -log(3), log(1), log(3), log(10), log(30), log(100), log(1000), log(10000)), yaxis.labels=c("1/30", "1/10", "1/3", "1", "3", "10", "30", "100", "1000", "10000"))

dev.off()



# Figure 2: unlimited SBF
pdf(file="Figure2.pdf", width=7, height=4, pointsize=10)
plot(sim, n.max=NA, xlim=c(10, 240), boundary=6, traj.selection="fixed", xextension=1.8)
dev.off()


# Figure 3: SBF+maxN
pdf(file="Figure3.pdf", width=7, height=4, pointsize=10)
plot(sim, n.max=60, xlim=c(10, 240), boundary=10, traj.selection="fixed", xextension=0.8)
dev.off()