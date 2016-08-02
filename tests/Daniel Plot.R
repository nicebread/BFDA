print(Sys.time())
library(reshape2)

sim0.7 <- BFDA.sim.ttest(d=0.7, n.min=10, n.max=250, boundary=Inf, stepsize=5, alternative="undirected", design="sequential", B=1000, verbose=TRUE, cores=2, r=sqrt(2)/2)

#load("../finalSims/BFDA.0.5.RData")

BFDA.analysis(sim0.7, boundary=6, n.max=200)

sim <- sim0.7$sim

BF.power <- sim %>% group_by(n) %>% summarise(
	BF3 = sum(logBF>log(3))/n(),
	BF6 = sum(logBF>log(6))/n(),
	BF10 = sum(logBF>log(10))/n(),
	BF30 = sum(logBF>log(30))/n(),
	BF100 = sum(logBF>log(100))/n()
)

BF.power.long <- melt(BF.power, id.vars="n")

ggplot(BF.power.long %>% filter(n<=150), aes(x=n, y=value, color=variable)) + geom_line() + ylab("Power = Prob( BF > cutoff | d = 0.7)") + xlab("Sample size per condition") + scale_color_discrete("Critical Bayes Factor") + ggtitle("Bayesian power for d = 0.7")

# Table
print(Sys.time())