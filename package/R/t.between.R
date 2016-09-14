## ======================================================================
## These are the simulation files for a two-group t-test
## ======================================================================

# ---------------------------------------------------------------------
# sample.function: return a data frame/ matrix with simulated raw data

# get two samples with a specific standardized mean differences
# ES = Cohen's d
sample.t.between <- function(n, ES, options.sample=NULL) {
	x <- rnorm(n, mean=0, sd=1)
	y <- rnorm(n, mean=0, sd=1) - ES
	return(cbind(x, y))
}

# ---------------------------------------------------------------------
# select.function: select a specified number of accumulating data from 
# the data frame/ matrix that was simulated with sample.function

select.t.between <- function(MAXSAMP, n) {
	return(MAXSAMP[1:n, ])
}


# ---------------------------------------------------------------------
# freq.test.function: return p.value, test statistic, and empirical ES

freq.test.t.between <- function(SAMP, alternative="directional") {

	if (alternative=="directional") {
		alt2 <- "greater"
	} else {
		alt2 <- "two.sided"
	}

	t1 <- t.test(SAMP[, 1], SAMP[, 2], var.equal=TRUE, alternative=alt2)

	return(list(
		statistic = t1$statistic,
		p.value = t1$p.value,
		emp.ES = 2*t1$statistic / sqrt(2*nrow(SAMP)-2)
	))
}

# ---------------------------------------------------------------------
# BF.test.function: return log(BF10)

BF.test.t.between <- function(SAMP, alternative="directional", freq.test=NULL, ...) {

	if (alternative=="directional") {
		nullInterval <- c(0, Inf)
	} else {
		nullInterval <- NULL
	}

	t.value <- 

	# suppress the "t is large; approximation invoked" message
	suppressMessages({
		t1 <- BayesFactor::ttest.tstat(freq.test$statistic, nrow(SAMP), nrow(SAMP), nullInterval=nullInterval, simple=TRUE, ...)
	})
	
	# returns the log(BF10)
	return(log(t1))
}