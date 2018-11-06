## ======================================================================
## These are the simulation files for a paired t-test
## ======================================================================

# ---------------------------------------------------------------------
# sample.function: return a data frame/ matrix / vector with simulated raw data

# get two samples with a specific standardized mean differences
# ES = dz (standardized difference scores)
sample.t.paired <- function(n, ES, options.sample=NULL) {
	rnorm(n, ES, sd=1)
}

# ---------------------------------------------------------------------
# select.function: select a specified number of accumulating data from 
# the data frame/ matrix that was simulated with sample.function

select.t.paired <- function(MAXSAMP, n) {
	return(MAXSAMP[1:n])
}


# ---------------------------------------------------------------------
# freq.test.function: return p.value, test statistic, and empirical ES

freq.test.t.paired <- function(SAMP, alternative=NULL, options.sample=NULL) {

	if (alternative=="directional") {
		alt2 <- "greater"
	} else {
		alt2 <- "two.sided"
	}

	t1 <- t.test(SAMP, mu=0, alternative=alt2)

	# see http://journal.frontiersin.org/article/10.3389/fpsyg.2013.00863/full

	# must returns these values
	return(list(
		statistic = t1$statistic,
		p.value = t1$p.value,
		emp.ES = t1$statistic / sqrt(length(SAMP))
	))
}

# ---------------------------------------------------------------------
# BF.test.function: return log(BF10)

BF.test.t.paired <- function(SAMP, alternative="directional", freq.test=NULL, ...) {

	if (alternative=="directional") {
		nullInterval <- c(0, Inf)
	} else {
		nullInterval <- NULL
	}

	# suppress the "t is large; approximation invoked" message
	suppressMessages({
		t1 <- BayesFactor::ttest.tstat(freq.test$statistic, n1=length(SAMP), n2=0, nullInterval=nullInterval, simple=TRUE, ...)
	})
	
	# returns the log(BF10)
	return(log(t1))
}