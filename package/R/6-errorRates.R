# TODO: This function is a stub!
errorRates <- function(BPA.H1, BPA.H0, n, criterion=c(1/3, 3), xlim=NA, noSplit=FALSE, cex=1.2, cex.axis=1) {

	# BPA.H1 <- BPA.0.5
	# BPA.H0 <- BPA.0
	# criterion = c(1/6, 6)
	
	# just in case: order criterion
	criterion <- sort(criterion)
	logCriterion <- log(criterion)

	n.max <- min(max(BPA.H1$sim$n), max(BPA.H0$sim$n))
	BPA.H1$sim <- BPA.H1$sim %>% filter(n <= n.max)
	BPA.H0$sim <- BPA.H0$sim %>% filter(n <= n.max)

	# ---------------------------------------------------------------------
	# Plot: x-axis = n, y-axis = error rate, fixed boundary
	
	# WE = weak evidence
	underH1 <- BPA.H1$sim %>% group_by(n) %>% summarize(
		FNE = sum(logBF < logCriterion[1])/n(),
		WE  = sum(logBF > logCriterion[1] & logBF < logCriterion[2])/n(),
		CE  = sum(logBF > logCriterion[2])/n()
	)
	
	underH0 <- BPA.H0$sim %>% group_by(n) %>% summarize(
		FPE = sum(logBF > logCriterion[2])/n(),
		WE  = sum(logBF > logCriterion[1] & logBF < logCriterion[2])/n(),
		CE  = sum(logBF < logCriterion[1])/n()
	)	

	p0 <- ggplot(underH0, aes(x=n, y=FPE)) + geom_smooth(span=.2, se=FALSE, color="darkred") + theme_bw()
	
	pWE <- ggplot(underH1, aes(x=n, y=WE)) + geom_smooth(span=.2, se=FALSE, color="darkgreen") + geom_smooth(data=underH0, aes(y=WE), span=.2, se=FALSE, color="darkred") + theme_bw()
	
	pCE <- ggplot(underH1, aes(x=n, y=WE)) + geom_smooth(span=.2, se=FALSE, color="darkgreen") + geom_smooth(data=underH0, aes(y=WE), span=.2, se=FALSE, color="darkred") + theme_bw()
	
	p1 <- ggplot(underH1, aes(x=n, y=FNE)) + geom_smooth(span=.2, se=FALSE, color="darkgreen") + theme_bw()
	
	library(gridExtra)
	grid.arrange(p0, p1, pWE)
}
	