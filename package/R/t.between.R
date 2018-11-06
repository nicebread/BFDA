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

freq.test.t.between <- function(SAMP, alternative=NULL) {

	t1 <- t.test(SAMP[, 1], SAMP[, 2], var.equal=TRUE, alternative=alternative)

	return(list(
		statistic = t1$statistic,
		p.value = t1$p.value,
		emp.ES = 2*t1$statistic / sqrt(2*nrow(SAMP)-2)
	))
}

# ---------------------------------------------------------------------
# Check definition of prior

prior.check.t.between <- function(prior=NULL){
  
  if(!is.list(prior)){
    if(!is.null(prior) == TRUE) {
      stop("Argument prior needs to be specified as a list.")
    } else {
      prior <- list("Cauchy", list(prior.location = 0, prior.scale = sqrt(2)/2))
    }}
  
  match.arg(prior[[1]], c("Cauchy", "t", "normal"))
  
  if(prior[[1]] %in% c("Cauchy", "t")){
    if(is.null(prior[[2]][["prior.location"]])){
      warning("Prior location not defined. Default specification will be used.")
      prior[[2]][["prior.location"]] <- 0
    } 
    if(is.null(prior[[2]][["prior.scale"]])){
      warning("Prior scale not defined. Default specification will be used.")
      prior[[2]][["prior.scale"]] <- sqrt(2)/2
    }
  }
  
  if(prior[[1]] == "Cauchy"){
    if(any(!is.element(names(prior[[2]]), c("prior.location", "prior.scale")))){
      warning("Cauchy distribution only takes parameters prior.location and prior.scale.")
    }
  } else if (prior[[1]] == "t"){
    if(any(!is.element(names(prior[[2]]), c("prior.location", "prior.scale", "prior.df")))){
      warning("t distribution only takes parameters prior.location, prior.scale, and prior.df.")
    }
    if(is.null(prior[[2]][["prior.df"]])) {
      warning("Prior degrees of freedom not defined. Default specifications will be used.")
      prior[[2]][["prior.df"]] <- 1
    }
  } else if (prior[[1]] == "normal") {
    if(any(!is.element(names(prior[[2]]), c("prior.mean", "prior.variance")))){
      warning("Normal distribution only takes parameters prior.mean and prior.variance.")
    }
    if(is.null(prior[[2]][["prior.mean"]])){
      warning("Prior mean not defined. Default specification will be used.")
      prior[[2]][["prior.mean"]] <- 0
    } 
    if(is.null(prior[[2]][["prior.variance"]])){
      prior[[2]][["prior.variance"]] <- 1
      warning("Prior variance not defined. Default specification will be used.")
    }
  }
  return(prior)
}


# ---------------------------------------------------------------------
# BF.test.function: return log(BF10)

BF.test.t.between <- function(SAMP, alternative=NULL, freq.test=NULL, prior=NULL, ...) {

	if (alternative=="greater") {
		alt2 <- "BFplus0"
	} else if (alternative=="two.sided"){
		alt2 <- "BF10"
	} else if (alternative=="less"){
	  alt2 <- "BFmin0"
	}
  
  
  if(prior[[1]] == "Cauchy"){
    
    t1 <- bf10_t(t = as.numeric(freq.test$statistic), n1 = nrow(SAMP), n2 = nrow(SAMP), independentSamples = TRUE, prior.location = prior[[2]][["prior.location"]],
           prior.scale = prior[[2]][["prior.scale"]], prior.df = 1)
    
  } else if (prior[[1]] == "t") {
    
    t1 <- bf10_t(t = as.numeric(freq.test$statistic), n1 = nrow(SAMP), n2 = nrow(SAMP), independentSamples = TRUE, prior.location = prior[[2]][["prior.location"]],
           prior.scale = prior[[2]][["prior.scale"]], prior.df = prior[[2]][["prior.df"]])
    
  } else if (prior[[1]] == "normal") {
  
      
    t1 <- bf10_normal(t = as.numeric(freq.test$statistic), n1=nrow(SAMP), n2 = nrow(SAMP), independentSamples = TRUE,
                      prior.mean = prior[[2]][["prior.mean"]], prior.variance = prior[[2]][["prior.variance"]])
  }
  
	# returns the log(BF10)
	return(as.numeric(log(t1[[alt2]])))
}
