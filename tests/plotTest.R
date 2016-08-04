## ======================================================================
## PLot the three plots for the paper
## ======================================================================

# ---------------------------------------------------------------------
#  Figure 1: fixed-n

load("../finalSims/BFDA.0.5.prior.RData")
load("../finalSims/BFDA.0.4.prior.RData")
load("../finalSims/BFDA.0.5.RData")
load("../finalSims/BFDA.0.RData")

(a <- BFDA.analyze(BFDA.0.5, boundary=6, n.max=50))
hist(a$endpoint.n)
quantile(a$endpoint.n, prob=.75)

BFDA.analyze(BFDA.0, boundary=6, n.max=50)


pdf(file="Figure1.pdf", width=7, height=4, pointsize=10)

par(oma = c(0, 0, 4, 0), mar=c(4, 2.3, 1, 1))

# create a layout with an empty column in the middle (relative width=1), 
# and two columns left and right (relative width=7)
l1 <- layout(matrix(1:6, ncol=3), widths = c(7, 1, 7), heights = c(4, 4), respect = FALSE)
##  Draw Plots
compDist(BFDA.0.5.prior, BFDA.0, n=20, boundary=log(c(1/3, 3)), xlim=c(log(1/10), log(90000)), noSplit=TRUE, cex=0.8)
plot.new();plot.new(); # skip empty column in the middle
compDist(BFDA.0.5.prior, BFDA.0, n=100, boundary=log(c(1/3, 3)), xlim=c(log(1/10), log(90000)), noSplit=TRUE, cex=0.8)

##  Create an overall title.
mtext("n = 20, d = 0.5", line=2, outer = TRUE, adj=0.2)
mtext("n = 100, d = 0.5", line=2, outer = TRUE, adj=0.8)
dev.off()


# Comparison of fixed and diffused prior

par(oma = c(0, 0, 4, 0), mar=c(4, 2.3, 1, 1))

# create a layout with an empty column in the middle (relative width=1), 
# and two columns left and right (relative width=7)
l1 <- layout(matrix(1:6, ncol=3), widths = c(7, 1, 7), heights = c(4, 4), respect = FALSE)
##  Draw Plots
compDist(BFDA.0.5, BFDA.0, n=20, boundary=log(c(1/3, 3)), xlim=c(log(1/10), log(90000)), noSplit=TRUE, cex=0.8)
plot.new();plot.new(); # skip empty column in the middle
compDist(BFDA.0.5.prior, BFDA.0, n=20, boundary=log(c(1/3, 3)), xlim=c(log(1/10), log(90000)), noSplit=TRUE, cex=0.8)

##  Create an overall title.
mtext("n = 40, d = 0.5", line=2, outer = TRUE, adj=0.2)
mtext("n = 40, d = N(0.5, SD=0.1)", line=2, outer = TRUE, adj=0.8)


BFDA.analyze(BFDA.0.4.prior, boundary=6, n.max=100)
SSD(BFDA.0.5, power=.90)
SSD(BFDA.0.5.prior, power=.90)

plot(BFDA.0.5, boundary=10, n.max=NA)



# ---------------------------------------------------------------------
#  Figure 2: open-ended SBF

pdf(file="Figure2.pdf", width=7, height=4, pointsize=10)
plot(BFDA.0.5.prior, n.max=NA, xlim=c(10, 180), boundary=6, xextension=2.2, dens.amplification=1.5, ylim=c(log(1/11), log(11)))
dev.off()


# Figure 3: SBF+maxN
pdf(file="Figure3.pdf", width=7, height=4, pointsize=10)
plot(BFDA.0.5.prior, n.max=60, boundary=6, traj.selection="fixed", xextension=1.5, cex.labels=1.1, cex.annotations=0.9)
dev.off()



power.NHST <- BFDA.0.5.prior$sim %>% group_by(n) %>% summarise(sig.prop = sum(p.value <= .05)/n())
print(power.NHST, n=Inf)

# ---------------------------------------------------------------------
# SSD plot test

# Correct sign of BF?
SSD(BFDA.0.5, power=.80, criterion=c(1/3, 3))
SSD(BFDA.0, alpha=.05, criterion=c(1/3, 3))


# symmetric boundaries at 1/3 and 3
SSD(sim, power=.90, criterion=c(1/3, 3))

# asymmetric boundaries at 0.5 and 10
SSD(sim, power=.80, criterion=c(0.5, 10))

# decrease Type-II errors with stronger
SSD(sim, power=.80, criterion=c(1/4, 10))


SSD(sim, power=.50, criterion=c(1/3, 3))

# Prepare lean ES=0 object
load("../finalSims/sim.0.RData")
BFDA0 <- sim2
BFDA0$sim <- BFDA0$sim[, c("d", "boundary", "n", "logBF", "d.emp", "id")]
