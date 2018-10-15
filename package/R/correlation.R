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

freq.test.correlation <- function(SAMP, alternative=NULL) {

	t1 <- cor.test(SAMP[, 1], SAMP[, 2], alternative=alternative)

	return(list(
		statistic = t1$statistic,
		p.value = t1$p.value,
		emp.ES = t1$estimate
	))
}

# ---------------------------------------------------------------------
# Check definition of prior

prior.check.correlation <- function(prior=NULL) {
  
  if(!is.list(prior)){
    if(!is.null(prior) == TRUE) {
      stop("Argument prior needs to be specified as a list.")
    } else {
      prior <- list("stretchedbeta", list(prior.alpha=1))
    }}
  
  match.arg(prior[[1]], "stretchedbeta")
  
  if(names(prior[[2]]) != "prior.alpha") {
    warning("Stretched beta prior only takes parameter prior.alpha. Default value will be applied.")
    prior[[2]][["prior.alpha"]] <- 1
  }
  
  if(prior[[2]][["prior.alpha"]] > 2 | prior[[2]][["prior.alpha"]] < 0 ) stop("Prior.alpha must be between 0 and 2 ")

  return(prior)
}

# ---------------------------------------------------------------------
# BF.test.function: return log(BF10)

BF.test.correlation <- function(SAMP, alternative=NULL, freq.test=NULL, prior=NULL, ...) {
  
  if(alternative == "greater" | alternative == "less") {
    alt2 <- "directional"
  } else if (alternative == "two.sided"){
    alt2 <- "two.sided"
  }

	if (alt2=="directional") {
		t1 <- corBF1sided(nrow(SAMP), freq.test$emp.ES, alpha=prior[[2]][["prior.alpha"]], ...)
	} else if (alt2 == "two.sided"){
		t1 <- corBF2sided(nrow(SAMP), freq.test$emp.ES, alpha=prior[[2]][["prior.alpha"]], ...)
	}
	
	# returns the log(BF10)
	return(as.numeric(log(t1)))
}

