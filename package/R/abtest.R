## =============================================================================
## These are the simulation files for the AB test
## =============================================================================

# ------------------------------------------------------------------------------
# sample.function: return a data frame/ matrix with simulated raw data

# findp2: abtest sample helper function to determine the second proportion depending on a given first proportion 

findp2 <- function(ES, p1, effecttype){
  
  if(effecttype == "OR"){
    num <- ES*(p1/(1-p1))
    denom <- 1+ES*(p1/(1-p1))
    p2 <- num/denom
  } else if(effecttype == "RR"){
    p2 <- p1*ES
  } else if(effecttype == "AR"){
    p2 <- p1+ES
  }
  
  return(p2)
}

#' sample.abtest
#' get two samples with a specified relation of underlying proportions
#' @importFrom stats rbinom
#' @param n sample size of the trial per condition
#' @param ES effect size (value of odds ratio, log odds ratio, absolute risk, relative risk)
#' @param options.sample list defining the type of the input effect size: \code{list(effecttype = ...)} with the options: \code{"OR"} (odds ratio), \code{"RR"} (relative risk), \code{"AR"} (absolute risk), \code{"logOR"} (log odds ratio)
#' @examples \dontrun{SAMP <- sample.abtest(n=100, ES=2, options.sample = list(effecttype = "OR"))}

sample.abtest <- function(n, ES, options.sample=NULL) {
  
  effecttype <- options.sample[["effecttype"]]
  
  if(effecttype == "OR" | effecttype == "logOR"){
    if(effecttype == "logOR") {ES <- exp(ES)} #replace log odds ratio by odds ratio
    p1 <- 0.5
    p2 <- findp2(ES = ES, p1 = p1, effecttype = "OR")
  }
  
  if(effecttype == "RR"){
    findp1RR <- function(RR){
      i <- 2
      p1 <- 0.5
      while(1/i * RR >= 1){
        p1 <- 1/(i+1)
        i <- i+1
      }
      return(p1)
    }
    p1 <- findp1RR(RR = ES)
    p2 <- findp2(ES = ES, p1 = p1, effecttype = "RR")
  }
  
  if(effecttype == "AR"){
    if(ES >= 1) stop("Absolute risk can only take values between 0 and 1")
    findp1AR <- function(AR){ # findp1: Find two values whose sum is smaller than one
      i <- 2
      p1 <- 0.5
      while(1/i + AR >= 1){
        p1 <- 1/(i+2)
        i <- i+1
        }
      return(p1)
    }
    p1 <- findp1AR(AR = ES)
    p2 <- findp2(ES = ES, p1 = p1, effecttype = "AR")
  }
  
  x <- stats::rbinom(n, size = 1, prob = p1)
  y <- stats::rbinom(n, size = 1, prob = p2)
  
  return(cbind(x, y))
}

# ------------------------------------------------------------------------------
# select.function: select a specified number of accumulating data from 
# the data frame/ matrix that was simulated with sample.function

select.abtest <- function(MAXSAMP, n) {
  return(MAXSAMP[1:n, ])
}

# ------------------------------------------------------------------------------
#' freq.test.function: return p.value, test statistic, and empirical ES
#' @param SAMP (reduced) sample
#' @param alternative direction of the alternative hypothesis
#' @param options.sample list(effecttype = ...)
#' @importFrom stats prop.test

freq.test.abtest <- function(SAMP, alternative=NULL, options.sample=NULL){
  
  t1 <- suppressWarnings(prop.test(x = as.numeric(colSums(SAMP)),
                         n = rep(nrow(SAMP), 2),
                         alternative = alternative)
                  )
  
  p1 <- t1$estimate[1]
  p2 <- t1$estimate[2]
  
  if(options.sample[["effecttype"]] == "OR"){
    emp.ES <- (p2/(1-p2))/(p1/(1-p1))
  } else if (options.sample[["effecttype"]] == "RR"){
    emp.ES <- p2/p1
  } else if (options.sample[["effecttype"]] == "AR"){
    emp.ES <- p2-p1
  } else if (options.sample[["effecttype"]] == "logOR"){
    emp.ES <- log((p2/(1-p2))/(p1/(1-p1)))
  }
  
  return(list(
    statistic = as.numeric(t1$statistic),
    p.value = as.numeric(t1$p.value),
    emp.ES = as.numeric(emp.ES)
  ))
  
}

# ------------------------------------------------------------------------------
# Check definition of prior

prior.check.abtest <- function(prior=NULL) {
  
  if(!is.list(prior)){
    if(!is.null(prior) == TRUE) {
      stop("Argument prior needs to be specified as a list.")
    } else {
      prior <- list("normal", list(prior.mean = 0, prior.variance = 1))
    }}
  
  match.arg(prior[[1]], "normal")
  
  if(any(!is.element(names(prior[[2]]), c("prior.mean", "prior.variance")))){
    warning("Normal distribution only takes parameters prior.mean and prior.variance.")
  }
  
  if(is.null(prior[[2]][["prior.mean"]])){
    warning("Prior mean not defined. Default specification will be used.")
    prior[[2]][["prior.mean"]] <- 0
  } 
  
  if(is.null(prior[[2]][["prior.variance"]])){
    prior[[2]][["prior.variance"]] <- sqrt(2)/2
    warning("Prior variance not defined. Default specification will be used.")
  }
  
  return(prior)
}

# ------------------------------------------------------------------------------
# BF.test.function: return log(BF10)

BF.test.abtest <- function(SAMP, alternative=NULL, freq.test=NULL, prior=NULL, ...) {
  
  if(alternative == "greater") {
    alt2 <- "bfplus0"
  } else if (alternative == "two.sided"){
    alt2 <- "bf10"
  } else if (alternative == "less"){
    alt2 <- "bfminus0"
  }
  
  testdat <- list(y1 = colSums(SAMP)[1],
                  y2 = colSums(SAMP)[2],
                  n1 = nrow(SAMP),
                  n2 = nrow(SAMP))
  
  testprior <- list(mu_psi = prior[[2]][["prior.mean"]],
                    sigma_psi = prior[[2]][["prior.variance"]],
                    mu_beta = 0,
                    sigma_beta = 1)
  
  t1 <- as.numeric(ab_test(data = testdat, prior_par = testprior)$bf[[alt2]])
  
  return(as.numeric(log(t1)))
  
}
