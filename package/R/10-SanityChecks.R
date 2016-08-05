#BFDA <- corr.H1

# ---------------------------------------------------------------------
#  At each fixed n, the success rate (via p values) should equal the calculated power
# This is only valid in the fixed-ES case!

BFDA.sanityCheck <- function(BFDA) {
	
	if (var(BFDA$sim$true.ES)>0) stop("This sanity check only works for BFDA simulations with a fixed effect size!")
	
	# ---------------------------------------------------------------------
	# Sanity check for H1 simulations with fixed ES	

	if (BFDA$sim$true.ES[1] > 0) {
		library(pwr)

		print("Comparing analytical and simulated frequentist power analysis ...")
		power.comparison <- data.frame()
		for (n in unique(BFDA$sim$n)) {
			sim <- BFDA$sim[BFDA$sim$n == n,]
	
			alt <- ifelse(BFDA$settings$alternative == "undirected", "two.sided", "greater")
	
			power.analytically <- switch(BFDA$settings$type,
				"t.between" = {pwr.t.test(d=sim$true.ES[1], n=n, alternative=alt, type="two.sample")$power},
				"t.within" = {pwr.t.test(d=sim$true.ES[1], n=n, alternative=alt, type="one.sample")$power},
				"correlation" = {pwr.r.test(r=sim$true.ES[1], n=n, alternative=alt)$power}
			)
			power.sim <- sum(sim$p.value<.05)/length(sim$p.value)
			power.comparison <- rbind(power.comparison, data.frame(power.analytically, power.sim))
		}

		plot(power.comparison, pch=20)
		abline(a=0, b=1, lty="dotted")
	}

	# ---------------------------------------------------------------------
	# Sanity check for H0 simulations
	cat("
	Under H0, simulated p-values should be uniformly distributed at each sample size
	--> Show p-curve at each n
	")


	if (BFDA$sim$true.ES[1] == 0) {
		set.seed(0xBEEF)
		ggplot(BFDA$sim %>% filter(n %in% c(min(n), sample(unique(n), 4), max(n))), aes(x=p.value)) + geom_histogram() + facet_wrap(~n)
	}
	
}