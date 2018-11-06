## ======================================================================
## These are the simulation files for correlations
## ======================================================================

# ---------------------------------------------------------------------
# helpers

# returns a data frame of two variables which correlate
getBiCop <- function(n, rho, mar.fun=rnorm, x = NULL, ...) {
	if (!is.null(x)) {X1 <- x} else {X1 <- mar.fun(n, ...)}
	if (!is.null(x) & length(x) != n) warning("Variable x does not have the same length as n!")
	
	C <- matrix(rho, nrow = 2, ncol = 2)
	diag(C) <- 1
	
	C <- chol(C)

	X2 <- mar.fun(n)
	X <- cbind(X1, X2)

	# induce correlation (does not change X1)
	df <- X %*% C

	return(df)
}


# ---------------------------------------------------------------------
# sample.function: return a data frame/ matrix with simulated raw data

# get two variables with a certain (population) correlation
# ES = correlation
sample.correlation <- function(n, ES, options.sample=NULL) {
	getBiCop(n, ES)
}

# ---------------------------------------------------------------------
# select.function: select a specified number of accumulating data from 
# the data frame/ matrix that was simulated with sample.function

select.correlation <- function(MAXSAMP, n) {
	return(MAXSAMP[1:n, ])
}


# ---------------------------------------------------------------------
# freq.test.function: return p.value, test statistic, and empirical ES

freq.test.correlation <- function(SAMP, alternative=NULL, options.sample=NULL) {

	if (alternative=="directional") {
		alt2 <- "greater"
	} else {
		alt2 <- "two.sided"
	}

	t1 <- cor.test(SAMP[, 1], SAMP[, 2], alternative=alt2)

	return(list(
		statistic = t1$statistic,
		p.value = t1$p.value,
		emp.ES = t1$estimate
	))
}

# ---------------------------------------------------------------------
# BF.test.function: return log(BF10)

BF.test.correlation <- function(SAMP, alternative="directional", freq.test=NULL, ...) {

	if (alternative=="directional") {
		t1 <- corBF1sided(nrow(SAMP), freq.test$emp.ES, ...)
	} else {
		t1 <- corBF2sided(nrow(SAMP), freq.test$emp.ES, ...)
	}
	
	# returns the log(BF10)
	return(log(t1))
}