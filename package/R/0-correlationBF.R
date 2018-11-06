# Correlation BF code
#
# The following script contains functions for the computation of a Bayesian
# correlation test which are implemented in JASP (see http://jasp-stats.org).
# The functions were retrieved from the JASP GitHub page on October 16, 2018
# and contain the scripts from the following two sub-pages:
# https://github.com/jasp-stats/jasp-desktop/blob/development/JASP-Engine/JASP/R/common.R
# https://github.com/jasp-stats/jasp-desktop/blob/development/JASP-Engine/JASP/R/correlationbayesian.R
#
# Copyright (C) 2013-2018 University of Amsterdam

# ------------------------------------------------------------------------------
# Helper function for Bayesian correlation

isTryError <- function(obj){
  if (is.list(obj)){
    return(any(sapply(obj, function(obj) {
      inherits(obj, "try-error")
    }))
    )
  } else {
    return(any(sapply(list(obj), function(obj){
      inherits(obj, "try-error")
    })))
  }
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Here starts the Bayesian correlation

# References for Bayesian correlation:
#
# Ly, A., Verhagen, A. J. & Wagenmakers, E.-J. (2016). Harold Jeffreys's Default
# Bayes Factor Hypothesis Tests: Explanation, Extension, and Application in Psychology.
# Journal of Mathematical Psychology, 72, 19-32.
# https://www.sciencedirect.com/science/article/pii/S0022249615000383
#
# Ly, A., Marsman, M., Wagenmakers, E.-J. (2018). Analytic Posteriors for 
# Pearson’s Correlation Coefficient. Statistica Neerlandica, 72(1), 4-13.
# https://onlinelibrary.wiley.com/doi/full/10.1111/stan.12111

## Help functions ------------------------------------------------------------
# 0.1 Prior specification Pearson's Rho
.excludePairwiseCorData <- function(v1, v2) {
  # To exclude the data pairwise
  #
  screenedData <- list(v1=v1, v2=v2)
  
  removeIndex1 <- which(is.na(v1))
  removeIndex2 <- which(is.na(v2))
  removeIndex <- unique(c(removeIndex1, removeIndex2))
  if (length(removeIndex) > 0) {
    screenedData$v1 <- v1[-(removeIndex)]
    screenedData$v2 <- v2[-(removeIndex)]
  }
  
  return(screenedData)
}

.stretchedBeta <- function(rho, alpha, beta) {
  result <- 1/2*dbeta((rho+1)/2, alpha, beta)
  return(result)
}


.priorRho <- function(rho, kappa=1) {
  .stretchedBeta(rho, 1/kappa, 1/kappa)
}

.priorRhoPlus <- function(rho, kappa=1) {
  nonNegativeIndex <- rho >=0
  lessThanOneIndex <- rho <=1
  valueIndex <- as.logical(nonNegativeIndex*lessThanOneIndex)
  result <- rho*0
  result[valueIndex] <- 2*.priorRho(rho[valueIndex], kappa)
  return(result)
}

.priorRhoMin <- function(rho, kappa=1) {
  negativeIndex <- rho <=0
  greaterThanMinOneIndex <- rho >= -1
  valueIndex <- as.logical(negativeIndex*greaterThanMinOneIndex)
  result <- rho*0
  result[valueIndex] <- 2*.priorRho(rho[valueIndex], kappa)
  return(result)
}

# 1.0. Built-up for likelihood functions
.aFunction <- function(n, r, rho) {
  hyperTerm <- Re(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=(1/2), z=(r*rho)^2))
  result <- (1-rho^2)^((n-1)/2)*hyperTerm
  return(result)
}

.bFunction <- function(n, r, rho) {
  hyperTerm <- Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=(3/2), z=(r*rho)^2))
  logTerm <- 2*(lgamma(n/2)-lgamma((n-1)/2))+((n-1)/2)*log(1-rho^2)
  result <- 2*r*rho*exp(logTerm)*hyperTerm
  return(result)
}

.hFunction <- function(n, r, rho) {
  result <- .aFunction(n=n, r=r, rho) + .bFunction(n=n, r=r, rho)
  return(result)
}

.hFunctionCombined <- function(nOri, rOri, nRep, rRep, rho) {
  result <- .hFunction(n=nOri, r=rOri, rho)*.hFunction(n=nRep, r=rRep, rho)
  return(result)
}

.hFunctionCombinedTwoSided <- function(nOri, rOri, nRep, rRep, rho) {
  result <- .aFunction(n=nOri, r=rOri, rho)*.aFunction(n=nRep, r=rRep, rho) +
    .bFunction(n=nOri, r=rOri, rho)*.bFunction(n=nRep, r=rRep, rho)
  return(result)
}

.hJeffreysApprox <- function(n, r, rho) {
  result <- ((1 - rho^(2))^(0.5*(n - 1)))/((1 - rho*r)^(n - 1 - 0.5))
  return(result)
}

# 1.1 Explicit marginal likelihood functions
.m0MarginalLikelihood <- function(s, t, n) {
  logTerm <- 2*lgamma(0.5*(n-1))
  result <- 1/4*n^(0.5*(1-2*n))*pi^(1-n)*(s*t)^(1-n)*exp(logTerm)
  return(result)
}

.m1MarginalLikelihoodNoRho <- function(s, t, n, r, rho) {
  return(.m0MarginalLikelihood(s, t, n)*
           (.aFunction(n=n, r=r, rho)+.bFunction(n=n, r=r, rho)))
}

#
# 2.1 Two-sided main Bayes factor ----------------------------------------------
.bf10Exact <- function(n, r, kappa=1) {
  # Ly et al 2015
  # This is the exact result with symmetric beta prior on rho
  # with parameter alpha. If kappa = 1 then uniform prior on rho
  if (n <= 2) {
    return(1)
  } else if (any(is.na(r))) {
    return(NA)
  }

  checkR <- abs(r) >= 1 # check whether |r| >= 1
  if (kappa >= 1 && n > 2 && checkR) {
    return(Inf)
  }
  #logHyperTerm <- log(hypergeo::hypergeo(((n-1)/2), ((n-1)/2), ((n+2/kappa)/2), r^2))
  logHyperTerm <- log(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=((n+2/kappa)/2), z=r^2))
  logResult <- log(2^(1-2/kappa))+0.5*log(pi)-lbeta(1/kappa, 1/kappa)+
    lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)+logHyperTerm
  realResult <- exp(Re(logResult))
  
  if (realResult < 0) {
    return(NA)
  }
  
  # Fail
  return(realResult)
}

# 2.2 Two-sided secondairy Bayes factor
.bf10JeffreysIntegrate <- function(n, r, kappa=1) {
  # Jeffreys' test for whether a correlation is zero or not
  # Jeffreys (1961), pp. 289-292
  # This is the exact result, see EJ
  ##
  if (n <= 2) {
    return(1)
  } else if ( any(is.na(r)) ) {
    return(NA)
  }
  
  # TODO: use which
  if (n > 2 && abs(r)==1) {
    return(Inf)
  }
  
  hyperTerm <- Re(hypergeo::genhypergeo(U=c((2*n-3)/4, (2*n-1)/4), L=(n+2/kappa)/2, z=r^2))
  logTerm <- lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)-lbeta(1/kappa, 1/kappa)
  result <- sqrt(pi)*2^(1-2/kappa)*exp(logTerm)*hyperTerm
  
  if (result < 0) {
    return(NA)
  }
  return(result)
}


# 2.3 Savage-Dickey beta approximation
.bfSavageDickeyBetaData <- function(n, r, kappa=1, rho0=0) {
  # Savage-Dickey based on a beta approximation
  #
  #
  fit <- .posteriorBetaParameters(n=n, r=r, kappa=kappa)
  return(.bfSavageDickeyBeta(a=fit$betaA, b=fit$betaB, kappa=kappa, rho0=rho0))
}

.bfSavageDickeyOneSidedAdapt <- function(bf10, a, b, kappa=1, rho0=0) {
  return(.bfSavageDickeyBeta(bf10=bf10, a=a, b=b, kappa))
}


.bfSavageDickeyBeta <- function(a, b, kappa=1, rho0=0, bf10=NULL) {
  # Savage-Dickey based on a beta approximation
  # Default failure for infinite bf10
  # Depending on bf10 define the result list
  
  if (is.null(bf10)) {
    result <- list(bf10=NA, bfPlus0=NA, bfMin0=NA, betaA=a, betaB=b)
    
    savageDickeyNumerator <- .priorRho(rho0, kappa)
    savageDickeyDenominator <- .stretchedBeta(rho0, alpha=a, beta=b)
    bf10 <- try(savageDickeyNumerator/savageDickeyDenominator)
    
    if (isTryError(bf10)) {
      # NAs
      return(result)
    } else {
      result$bf10 <- bf10
    }
  } else {
    result <- list(bf10=bf10, bfPlus0=NA, bfMin0=NA, betaA=a, betaB=b)
  }
  
  if (is.na(bf10)) {
    # Failure
    return(result)
  }
  
  if (bf10 < 0) {
    # Total failure it's true
    return(result)
  }
  
  if (is.infinite(bf10)) {
    # .bfCorrieKernel gives the default values for the one sided ones
    #
    result$bf10 <- Inf
    return(result)
  }
  
  if (is.finite(bf10)) {
    # bf10 is finite, now calculate one-sided stuff
    #
    
    result$bf10 <- bf10
    result$twoSidedTooPeaked <- FALSE
    
    leftProportion <- stats::pbeta(1/2, shape1=a, shape2=b)
    
    if (is.na(leftProportion)) {
      result <- utils::modifyList(result, list(bfPlus0=NA, bfMin0=NA, minSidedTooPeaked=TRUE, plusSidedTooPeaked=TRUE))
      return(result)
    }
    
    if (leftProportion > 0 && leftProportion < 1) {
      result$plusSidedTooPeaked <- FALSE
      result$minSidedTooPeaked <- FALSE
      
      result$bfMin0 <- 2*bf10*leftProportion
      result$bfPlus0 <- 2*bf10*(1-leftProportion)
    } else {
      rightProportion <- stats::pbeta(1/2, shape1=a, shape2=b, lower.tail=FALSE)
      
      if (!is.na(rightProportion) && rightProportion < 1) {
        result$plusSidedTooPeaked <- FALSE
        result$minSidedTooPeaked <- FALSE
        
        result$bfMin0 <- 2*bf10*(1-rightProportion)
        result$bfPlus0 <- 2*bf10*rightProportion
      } else if (leftProportion >= 1) {
        result$bfMin0 <- 2*bf10
        result$bfPlus0 <- 0
      } else if (rightProportion >= 1) {
        result$bfMin0 <- 0
        result$bfPlus0 <- 2*bf10
      } else {
        # TODO: either leftProportion <= 0 or rightProportion < = 0
        # this shouldn't happen, but actually I have no idea
        return(result)
      }
    }
  }
  return(result)
}

# 2.4 The Marsman MH sampler


.logTarget <- function(rho, n, r, kappa=1) {
  
  (1/kappa+(n-1)/2)*log(1-rho^2)-(2*n-3)/2*log(1-rho*r)
}


.logProposal <- function(rho, n, r, z=NULL) {
  # The proposal is defined on Fisher hyperbolic tangent transform of rho, 
  # the Jacobian is absorbed in .logTarget
  if (!is.null(z)) {
    return(-(n-3)/2*(z-atanh(r))^2)
  } else {
    -(n-3)/2*(atanh(rho)-atanh(r))^2
  }
}

.marsmanMHSampler <- function(n, r, kappa=1, nIters=50000) {
  rhoMetropolisChain <- numeric(nIters)
  
  if (n <= 3) {
    std <- 1
  } else {
    std <- 1 / sqrt(n - 3)
  }
  
  zCandidates <- rnorm(nIters, mean=atanh(r), sd=std)
  rhoCandidates <- tanh(zCandidates)
  logTargetCandidates <- .logTarget(rho=rhoCandidates, n=n, r=r, kappa=kappa)
  logPropCandidates <- .logProposal(z=zCandidates, n=n, r=r)
  acceptMechanism <- runif(nIters)
  candidateAcceptance <- numeric(nIters)
  
  rhoCurrent <- r
  
  for (iter in 1:nIters) {
    zCurrent <- atanh(rhoCurrent)
    
    candidateAcceptance[iter] <- logTargetCandidates[iter]+.logProposal(z=zCurrent, n=n, r=r)-
      (.logTarget(rho=rhoCurrent, n=n, r=r, kappa=kappa)+logPropCandidates[iter])
    
    if (log(acceptMechanism[iter]) <= candidateAcceptance[iter]) {
      # Accept candidate and update rhoCurrent for next iteration
      rhoCurrent <- rhoCandidates[iter]
    } 
    rhoMetropolisChain[iter] <- rhoCurrent
  }
  
  acceptanceRate <- length(unique(rhoMetropolisChain))/nIters
  
  metropolisVar <- var(rhoMetropolisChain)/2^2
  metropolisMean <- mean((rhoMetropolisChain+1)/2)
  
  mhFit <- .betaParameterEstimates(metropolisMean, metropolisVar)
  mhFit[["acceptanceRate"]] <- acceptanceRate
  return(mhFit)
}

# 2.5. Two-sided asymptotic approximation of the Bayes factor
.bf10JeffreysApprox <- function(n, r) {
  #Jeffreys' test for whether a correlation is zero or not
  #Jeffreys (1961), pp. 291 Eq. 14
  #
  result <- list(bf10=NA, bfPlus0=NA, bfMin0=NA)
  
  if (n <= 2) {
    return(1)
  } else if ( any(is.na(r)) ) {
    return(NA)
  }
  # TODO: use which
  if (n > 2 && abs(r)==1) {
    return(Inf)
  }
  
  bf01 <- ((2*n-3)/pi)^(.5)*(1-r^2)^((n-4)/2)
  bf10 <- 1/bf01
  
  if (bf10 < 0) {
    return(result)
  }
  
  if (is.na(bf10)) {
    return(result)
  }
  
  if (is.finite(bf10)) {
    result$bf10 <- bf10
    result$twoSidedTooPeaked <- FALSE
    
    return(result)
  }
  
  if (is.infinite(bf10)) {
    # Note: Check extreme
    if (r >= 0) {
      result$bfPlus0 <- Inf
      result$bfMin0 <- 0
    } else if (r < 0) {
      result$bfPlus0 <- 0
      result$bfMin0 <- Inf
    }
    
    # Posterior too peaked
    result$twoSidedTooPeaked <- TRUE
    result$plusSidedTooPeaked <- TRUE
    result$minSidedTooPeaked <- TRUE
    return(result)
  }
}

# 3.0 One-sided preparation ----------------------------------------------------
# For .bfPlus0Exact
# For .bfPlus0Exact
.mPlusExact <- function(n, r, kappa=1) {
  
  hyperTerm <- Re(hypergeo::genhypergeo(U=c(1, n/2, n/2), L=c(3/2, (2+kappa*(n+1))/(2*kappa)), z=r^2))
  logTerm <- 2*(lgamma(n/2)-lgamma((n-1)/2))-lbeta(1/kappa, 1/kappa)
  result <- 2^((3*kappa-2)/kappa)*kappa*r/(2+(n-1)*kappa)*exp(logTerm)*hyperTerm
  return(result)
}

# For .bfPlus0EJeffreysIntegrate
.mPlusJeffreysIntegrate <- function(n, r, kappa=1) {
  # Ly et al 2015
  # This is the exact result with symmetric beta prior on rho
  # This is the contribution of one-sided test
  #
  #
  hyperTerm <- Re(hypergeo::genhypergeo(U=c(1, (2*n-1)/4, (2*n+1)/4),
                                        L=c(3/2, (n+1+2/kappa)/2), z=r^2))
  logTerm <- -lbeta(1/kappa, 1/kappa)
  result <- 2^(1-2/kappa)*r*(2*n-3)/(n+2/kappa-1)*exp(logTerm)*hyperTerm
  return(result)
}



## Suit:
.bfHypergeo <- function(n, r, kappa=1, methodNumber=1, hyperGeoOverFlowThreshold=24) {
  # Outputs:
  #   list of bfs and the beta fits based on the exact form of the reduced likelihood,
  #   see Ly, Marsman and Wagenmakers (2017) "Analytic Posteriors for Pearson’s Correlation Coefficient".
  #
  # This is the exact result with symmetric beta prior on rho
  # with parameter alpha. If kappa = 1 then uniform prior on rho
  #
  
  result <- list(bf10=NA, bfPlus0=NA, bfMin0=NA)
  tempList <- .posteriorBetaParameters(n=n, r=r, kappa=kappa)
  result <- utils::modifyList(result, tempList)
  
  bf10 <- switch(methodNumber,
                 try(silent=TRUE, .bf10Exact(n=n, r=r, kappa=kappa)),
                 try(silent=TRUE, .bf10JeffreysIntegrate(n=n, r=r, kappa))
  )
  
  if (isTryError(bf10)) {
    # all NAs
    return(result)
  }
  
  if (is.na(bf10)) {
    # all NAs
    return(result)
  }
  
  if (bf10 <0) {
    # Total failure it's true
    return(result)
  }
  
  if (is.finite(bf10)) {
    # Store
    result$bf10 <- bf10
    result$twoSidedTooPeaked <- FALSE
    
    if (log(bf10) < hyperGeoOverFlowThreshold) {
      # No overflow, can use exact result
      switch(methodNumber,
             {
               bfPlus0 <- try(silent=TRUE, bf10 + .mPlusExact(n=n, r=r, kappa))
               bfMin0 <- try(silent=TRUE, bf10 + .mPlusExact(n, -r, kappa))
             },
             {
               bfPlus0 <- try(silent=TRUE, bf10 + .mPlusJeffreysIntegrate(n=n, r=r, kappa=kappa))
               bfMin0 <- try(silent=TRUE, bf10 + .mPlusJeffreysIntegrate(n=n, r=-r, kappa=kappa))
             }
      )
      
      if (is.finite(bfPlus0) && is.finite(bfMin0)) {
        tempList <- list(bfPlus0=bfPlus0, bfMin0=bfMin0, plusSidedTooPeaked=FALSE, minSidedTooPeaked=FALSE)
        result <- utils::modifyList(result, tempList)
        return(result)
      }
    }
    
    # bf10 either overflow, or exact one-sided bfs failed:
    # but bf10 is finite, so try Savage-Dickey adaptation for one-sided bfs
    tempList <- .bfSavageDickeyOneSidedAdapt(bf10, a=result$betaA, b=result$betaB, kappa=kappa)
    result <- utils::modifyList(result, tempList)
    return(result)
  }
  
  if (is.infinite(bf10)) {
    # Note: Information consistency check is done at a higher level
    #  For information consistent data it actually outputs the NA list
    #
    # Exact bf10 is infinite, perhaps due to hypergeometric, try numerical integrate
    twoSidedIntegrand <- switch(methodNumber, {
      function(x){.hFunction(n=n, r=r, x)*.priorRho(x, kappa=kappa)} }, {
        function(x){.hJeffreysApprox(n=n, r=r, x)*.priorRho(x, kappa=kappa)}
      })
    
    
    bf10 <- try(silent=TRUE, expr=integrate(twoSidedIntegrand, -1, 1)$value)
    result$bf10 <- bf10
    
    if (is.infinite(result$bf10)) {
      # .bfCorrieKernel gives the default values for the one sided ones
      #
      return(result)
    }
    
    if (is.finite(bf10)) {
      # Numerical integrated bf10 is finite
      result$twoSidedTooPeaked <- FALSE
      
      if (methodNumber==1) {
        plusSidedIntegrand <- function(x){.hFunction(n=n, r=r, x)*.priorRhoPlus(x, kappa=kappa)}
        minSidedIntegrand <- function(x){.hFunction(n=n, r=r, x)*.priorRhoMin(x, kappa=kappa)}
      }
      
      if (methodNumber==2) {
        plusSidedIntegrand <- function(x){.hJeffreysApprox(n=n, r=r, x)*.priorRhoPlus(x, kappa=kappa)}
        minSidedIntegrand <- function(x){.hJeffreysApprox(n=n, r=r, x)*.priorRhoMin(x, kappa=kappa)}
      }
      
      bfPlus0 <- try(silent=TRUE, expr=integrate(plusSidedIntegrand, 0, 1)$value)
      bfMin0 <- try(silent=TRUE, expr=integrate(minSidedIntegrand, -1, 0)$value)
      
      if (is.finite(bfPlus0) && is.finite(bfMin0)) {
        tempList <- list(bf10=bf10, bfMin0=bfMin0, bfPlus0=bfPlus0, plusSidedTooPeaked=FALSE, minSidedTooPeaked=FALSE)
        result <- utils::modifyList(result, tempList)
        return(result)
      }
      
      # Numerical procedure failed, but bf10 is finite, try savage-dickey adapt
      tempList <- .bfSavageDickeyOneSidedAdapt(bf10, a=result$betaA, b=result$betaB, kappa=kappa)
      result <- utils::modifyList(result, tempList)
    }
    
    # numerical bf10 is infinite
    #
    
    
    # Numerical bf10 is not finite nor infinite
    # All NAs
    return(result)
  }
  # Exact bf10 is not finite, nor infinite, nor NA, nor try-error No clue what happened here
  # All NAs
  return(result)
}

.bfCorrieKernel <- function(n, r, kappa=1, method="exact", ciValue=0.95, hyperGeoOverFlowThreshold=24) {
  # The idea is incremental when it comes to method numbers, if 1 doesn't work then go down.
  # In particular, when
  #	methodNumber=1: exact result Ly et al (2015)
  #	methodNumber=2: semi-exact result, based on approximation of the likelihood JeffreysExact, see Wagenmakers et al (2015) bathing
  #	methodNumber=3: Savage Dickey beta approximation
  #	methodNumber=4: Marsman's IMH sampler and report a posterior beta fit summarised by betaA, betaB
  #   methodNumber=5: Jeffreys asymptotic approximation of the BF
  #   methodNumber=6: invalid data, or point prior
  #
  #   hyperGeoOverFlowThreshold=25 implies that if log(bf10) > 24 that we use Savage-Dickey adaptation
  #   for the one-sided bfs.
  #
  # Ex. 	if we want confidence intervals and the bfs are based on method 1,
  #		then we can find betaA, betaB based on .posteriorBetaParameters
  # 		if this doens't work we use the Marsman sampler
  # Ex.	if we retrieve bfs from methodNumber 4, we know that we already have betaA and betaB,
  #		so we do not need to try .posteriorMean and .posteriorVariance
  # Ex.	if we retrieve bfs from methodNumber 6, we don't know anything about the posterior, so no report at all
  #
  
  tempList <- list()
  result <- list(n=n, r=r, kappa=kappa, bf10=NULL, bfPlus0=NULL, bfMin0=NULL, methodNumber=NULL, betaA=NULL, betaB=NULL,
                 twoSidedTooPeaked=NULL, plusSidedTooPeaked=NULL, minSidedTooPeaked=NULL,
                 ci=tempList, ciValue=ciValue, acceptanceRate=1)
  
  # When bf10 is na, modify result with naList
  #
  tooPeakedList <- list(twoSidedTooPeaked=TRUE, plusSidedTooPeaked=TRUE, minSidedTooPeaked=TRUE)
  
  naList <- append(list(bf10=NA, bfPlus0=NA, bfMin0=NA), tooPeakedList)
  
  # When the prior is trivial (null is alternative) or when the data is predictively matched
  #
  predictiveMatchingList <- list(bf10=1, bfPlus0=1, bfMin0=1, twoSidedTooPeaked=FALSE, plusSidedTooPeaked=FALSE, minSidedTooPeaked=FALSE, methodNumber=0)
  
  # Information consistent result
  #
  infoConsistentList <- append(list(bf10=Inf, methodNumber=0), tooPeakedList)
  
  # Note: If log(bf10) then use beta fit for the one sided bfs
  # The 24 is based on n=30, r=-0.3500786
  # hyperGeoOverFlowThreshold <- 24
  
  # Note: Data check
  #
  if (any(is.na(r)) ) {
    result$methodNumber <- 6
    result <- utils::modifyList(result, naList)
    return(result)
  }
  
  # Note: Data: OK
  # "No" prior, alternative model is the same as the null model
  # TODO: however this bound of 0.002 is arbitrarily chosen. I should choose this based on a trade off
  # between r and n, but it doesn't really matter.
  if (kappa <= 0.002) {
    result <- utils::modifyList(result, predictiveMatchingList)
    return(result)
  }
  
  checkR <- abs(r) >= 1 # check whether |r| >= 1
  if (n <= 2 || kappa==0) {
    result <- utils::modifyList(result, predictiveMatchingList)
    return(result)
  } else if (kappa >= 1 && n > 2 && checkR) {
    if (r > 0) {
      result$bfPlus0 <- Inf
      result$bfMin0 <- 0
    } else if (r <= 0) {
      result$bfPlus0 <- 0
      result$bfMin0 <- Inf
    }
    result <- utils::modifyList(result, infoConsistentList)
    return(result)
  }
  
  # Note: Define different methods and method number
  #
  if (method=="exact" || method==1) {
    result$methodNumber <- 1
    tempList <- .bfHypergeo(n=n, r=r, kappa=kappa, methodNumber=1, hyperGeoOverFlowThreshold=hyperGeoOverFlowThreshold)
    result <- utils::modifyList(result, tempList)
  }
  
  if (method=="jeffreysIntegrate" || method==2) {
    result$methodNumber <- 2
    tempList <- .bfHypergeo(n=n, r=r, kappa=kappa, methodNumber=2, hyperGeoOverFlowThreshold=hyperGeoOverFlowThreshold)
    result <- utils::modifyList(result, tempList)
  }
  
  if (method=="savageDickeyBeta" || method==3) {
    result$methodNumber <- 3
    tempList <- .bfSavageDickeyBetaData(n=n, r=r, kappa=kappa)
    result <- utils::modifyList(result, tempList)
  }
  
  if  (method=="metropolisHastings" || method==4) {
    # We use the Marsman Sampler (c) here based on posterior model fit.
    result$methodNumber <- 4
    marsmanResult <- .marsmanMHSampler(n=n, r=r, kappa=kappa)
    result$acceptanceRate <- marsmanResult$acceptanceRate
    
    tempResult <- .bfSavageDickeyBeta(a=marsmanResult$betaA, b=marsmanResult$betaB, kappa=kappa)
    result <- utils::modifyList(result, tempResult)
  }
  
  if (method=="jeffreysApprox" || method==5) {
    result$methodNumber <- 5
    tempResult <- .bf10JeffreysApprox(n=n, r=r)
    result <- utils::modifyList(result, tempResult)
    return(result)
  }
  
  # Note: bf10: CHECK
  if (is.na(result$bf10)) {
    # Posterior not interesting
    result <- utils::modifyList(result, naList)
    return(result)
  }
  
  if (is.infinite(result$bf10)) {
    if (r >= 0) {
      result$bfPlus0 <- Inf
      result$bfMin0 <- 0
    } else if (r < 0) {
      result$bfPlus0 <- 0
      result$bfMin0 <- Inf
    }
    
    result <- utils::modifyList(result, tooPeakedList)
    return(result)
  }
  
  # Note: Calculate credible intervals
  #
  if (!is.null(ciValue)) {
    # Note: ciValue=NULL, speeds up the calculations for sequential analysis
    result$ci <- .computePearsonCredibleInterval(alpha=result$betaA, beta=result$betaB, ciValue=result$ciValue)
  }
  
  # Note: bfPlus0, bfMin0: CHECK
  #
  if (any(is.na(result$bfPlus0), is.na(result$bfMin0))) {
    # Note: bfPlus0, bfMin0: NOT ok
    # 	if one is NA, then both are NA
    result$bfPlus0 <- NA
    result$bfMin0 <- NA
    
    # Posterior not interesting
    result$plusSidedTooPeaked <- TRUE
    result$minSidedTooPeaked <- TRUE
    return(result)
  }
  
  if (any(c(result$bfPlus0, result$bfMin0)==0)) {
    # Note: bfPlus0, bfMin0: EXTREME
    # 	if one is extreme, so is the other
    if (result$bfPlus0==0) {
      result$bfPlus0 <- 0
      result$bfMin0 <- Inf
    } else if (result$bfMin0==0) {
      result$bfPlus0 <- Inf
      result$bfMin0 <- 0
    }
    
    # Posterior too peaked
    result$plusSidedTooPeaked <- TRUE
    result$minSidedTooPeaked <- TRUE
    return(result)
  }
  
  
  # Note: bfPlus0, bfMin0: CHECK COHERENCE:
  if (result$bfPlus0 > 1 && result$bfMin0 > 1 || any(c(result$bfPlus0, result$bfMin0)<0)) {
    if (r > 0) {
      # Note: Data: OK,
      # 		bf10: OK.
      #		bfPlus0: OK
      #		bfMin0: NOT ok
      #
      # bfMin0 is bigger than one due to overflow: bfMin0 = 2*bf10 - bfPlus0.
      # Example: 2*1.2.... 10^ 24 - 2.... 10^24 = 1... 10^12 (due to round off)
      #
      result$bfMin0 <- 10^(-317)
      result$bfPlus0 <- 2*result$bf10 - result$bfMin0
    } else if (r < 0) {
      # Note: Data: OK,
      # 		bf10: OK.
      #		bfPlus0: NOT ok
      #		bfMin0: OK
      result$bfPlus0 <- 10^(-317)
      result$bfMin0 <- 2*result$bf10 - result$bfPlus0
    }
  }
  return(result)
}

bfPearsonCorrelation <- function(n, r, kappa=1, ciValue=0.95, hyperGeoOverFlowThreshold=24) {
  # Wrapper around .bfCorrieKernel
  #
  result <- list(bf10=NA, bfPlus0=NA, bfMin0=NA)
  methodNumber <- 1
  
  while (any(is.na(c(result$bf10, result$bfPlus0, result$bfMin0)),
             is.infinite(c(result$bf10, result$bfPlus0, result$bfMin0)))
         && methodNumber <= 4) {
    # Note: Try all normal methods
    # 1. Exact
    # 2. semi-exact result
    # 3. Savage-Dickey beta approximation
    # 4. Marsman sampler
    
    result <- .bfCorrieKernel(n=n, r=r, kappa=kappa, method=methodNumber, ciValue=ciValue, hyperGeoOverFlowThreshold=hyperGeoOverFlowThreshold)
    #}
    methodNumber <- methodNumber+1
  }
  
  result$call <- paste0(".bfPearsonCorrelation(n=", n, ", r=", r, ", kappa=", kappa, ", ciValue=", ciValue, ", hyperGeoOverFlowThreshold=", hyperGeoOverFlowThreshold, ")")
  result$stat <- r
  return(result)
}


# 4.0 Posteriors to graph TODO: we have to think about this, different
# results, thus, also switching of the illustrations?
#


# 4.1 Two-sided
.posteriorRho <- function(n, r, rho, kappa=1) {
  if (!is.na(r) && !r==0) {
    return(1/.bf10Exact(n=n, r=r, kappa)*.hFunction(n=n, r=r, rho)*.priorRho(rho, kappa))
  } else if (!is.na(r) && r==0) {
    return(1/.bf10JeffreysIntegrate(n=n, r=r, kappa)*.hJeffreysApprox(n=n, r=r, rho)*.priorRho(rho, kappa))
  }
}

.posteriorRhoPlus <- function(n, r, rho, kappa=1) {
  if (!is.na(r) && !r==0) {
    return(1/.bfCorrieKernel(n=n, r=r, kappa, method="exact")$bfPlus0*.hFunction(n=n, r=r, rho)*.priorRhoPlus(rho, kappa))
  } else if (!is.na(r) && r==0) {
    return(1/.bfCorrieKernel(n=n, r=r, kappa, method="jeffreysIntegrate")$bfPlus0*.hJeffreysApprox(n=n, r=r, rho)*.priorRhoPlus(rho, kappa))
  }
}

.posteriorRhoMin <- function(n, r, rho, kappa=1) {
  if (!is.na(r) && !r==0) {
    return(1/.bfCorrieKernel(n=n, r=r, kappa, method="exact")$bfMin0*.hFunction(n=n, r=r, rho)*.priorRhoMin(rho, kappa))
  } else if (!is.na(r) && r==0) {
    return(1/.bfCorrieKernel(n=n, r=r, kappa, method="jeffreysIntegrate")$bfMin0*.hJeffreysApprox(n=n, r=r, rho)*.priorRhoMin(rho, kappa))
  }
  
}


.approximatePosteriorRho <- function(rho, n, r) {
  1/(1-rho^2)*stats::dnorm(atanh(rho), mean=atanh(r), sd=1/sqrt(n))
}

.approximatePosteriorRhoPlus <- function(rho, n, r) {
  (.approximatePosteriorRho(rho,n,r) * (rho>0)) / (stats::pnorm(0,mean=atanh(r),sd=1/sqrt(n),lower.tail = FALSE))
}

.approximatePosteriorRhoMin <- function(rho, n, r) {
  (.approximatePosteriorRho(rho,n,r) * (rho<0)) / (stats::pnorm(0,mean=atanh(r),sd=1/sqrt(n)))
}



# 4.2
.posteriorMean <- function(n, r, kappa=1) {
  
  
  if (n <= 2) {
    return(NA)
  } else if (any(is.na(r))) {
    return(NA)
  }
  # TODO: use which
  checkR <- abs(r) >= 1 # check whether |r| >= 1
  
  if (kappa >= 1 && n > 2 && checkR) {
    return(r)
  }
  #logHyperTerm <- log(hypergeo::hypergeo(((n-1)/2), ((n-1)/2), ((n+2/kappa)/2), r^2))
  hyperTerm1 <- Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=c((2+(n+2)*kappa)/(2*kappa)), z=r^2))
  hyperTerm2 <- Re(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=c((2+n*kappa)/(2*kappa)), z=r^2))
  
  logResult <- 2*(lgamma(n/2)-lgamma((n-1)/2))
  result <- (2*kappa*r)/(2+n*kappa)*exp(logResult)*hyperTerm1/hyperTerm2
  
  if (is.na(result) || abs(result) > 1) {
    return(r)
  } else {
    return(result)
  }
}


.posteriorSecondMoment <- function(n, r, kappa=1) {
  #
  #
  if (n <= 2) {
    return(NA)
  } else if (any(is.na(r))) {
    return(NA)
  }
  # TODO: use which
  checkR <- abs(r) >= 1 # check whether |r| >= 1
  
  if (kappa >= 1 && n > 2 && checkR) {
    return(r)
  }
  #log.hyper.term <- log(hypergeo::hypergeo(((n-1)/2), ((n-1)/2), ((n+2/kappa)/2), r^2))
  hyperTerm1 <- Re(hypergeo::genhypergeo(U=c(3/2, (n-1)/2, (n-1)/2),
                                         L=c(1/2, (2+(n+2)*kappa)/(2*kappa)), z=r^2))
  hyperTerm2 <- Re(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=c((2+n*kappa)/(2*kappa)), z=r^2))
  
  result <- kappa/(n*kappa+2)*hyperTerm1/hyperTerm2
  
  if (is.na(result) || result <= 0) {
    return(NA)
  } else {
    return(result)
  }
}

.posteriorVariance <- function(n, r, kappa=1) {
  
  
  result <- .posteriorSecondMoment(n,r,kappa)-(.posteriorMean(n,r,kappa))^2
  
  if (is.na(result) | result <= 0) {
    return(NA)
  } else {
    return(result)
  }
}

.betaParameterEstimates <- function(someMean, someVar) {
  # someMean \in (0, 1)
  # TODO: think about someMean = 0
  some.a <- someMean*(someMean*(1-someMean)/someVar-1)
  some.b <- (1-someMean)*(someMean*(1-someMean)/someVar-1)
  
  result <- list(betaA=some.a, betaB=some.b)
  return(result)
}

.posteriorBetaParameters <- function(n, r, kappa=1) {
  some.mu <- try((.posteriorMean(n=n, r=r, kappa)+1)/2)
  some.var <- try(.posteriorVariance(n=n, r=r, kappa)/2^2)
  
  if (is(some.mu, "try-error") || is(some.var, "try-error") || is.na(some.mu) || is.na(some.var)) {
    # TODO: Before doing this try the MH sampler
    return(list(betaA=NA, betaB=NA))
  } else {
    return(.betaParameterEstimates(some.mu, some.var))
  }
}

.computePearsonCredibleInterval <- function(alpha, beta, ciValue) {
  # Compute Pearson's correlation credible interval based on a beta fit
  #
  result <- list(twoSided=NA, minSided=NA, plusSided=NA)
  
  if (is.null(ciValue)) {
    return(result)
  }
  
  if (ciValue <= 0 || ciValue >= 1) {
    # Note: Don't modify ciValue <- solve this at a higher level
    return(result)
  }
  
  typeOne <- 1-ciValue
  excessLevel <- typeOne/2
  
  if (any(is.na(c(alpha, beta)), is.infinite(c(alpha, beta)))) {
    return(result)
  } else {
    # Note: Zero one refers to the problem on the (0, 1) rather than on (-1, 1)
    lowerCIZeroOne <- try(qbeta(excessLevel, alpha, beta))
    medianCIZeroOne <- try(qbeta(1/2, alpha, beta))
    upperCIZeroOne <- try(qbeta(1-excessLevel, alpha, beta))
    
    if (isTryError(list(lowerCIZeroOne, medianCIZeroOne, upperCIZeroOne))) {
      return(result)
    } else {
      # Note: This is simply an application of the definition of the stretched beta
      lowerCI <- 2*lowerCIZeroOne-1
      medianCI <- 2*medianCIZeroOne-1
      upperCI <- 2*upperCIZeroOne-1
    }
  }
  
  
  result$twoSided <- c(lowerCI, medianCI, upperCI)
  
  # One sided:
  result$minSided <- .computePearsonMinSidedCredibleInterval(alpha, beta, ciValue)
  # The problem is symmetric
  temp  <- .computePearsonMinSidedCredibleInterval(beta, alpha, ciValue)
  result$plusSided <- c(-temp[3], -temp[2], -temp[1])
  
  return(result)
}



.computePearsonMinSidedCredibleInterval <- function(alpha, beta, ciValue) {
  # Compute min sided Pearson's correlation credible interval based on a beta fit
  #
  result <- NA
  typeOne <- 1-ciValue
  excessLevel <- typeOne/2
  
  if (any(is.na(c(alpha, beta)), is.infinite(c(alpha, beta)))) {
    return(result)
  } else {
    leftArea <- pbeta(1/2, alpha, beta)
    lowerCIZeroOne <- try(qbeta(excessLevel*leftArea, alpha, beta))
    medianCIZeroOne <- try(qbeta(leftArea/2, alpha, beta))
    upperCIZeroOne <- try(qbeta((1-excessLevel)*leftArea, alpha, beta))
    
    if (isTryError(list(lowerCIZeroOne, medianCIZeroOne, upperCIZeroOne))) {
      return(result)
    } else {
      lowerCI <- 2*lowerCIZeroOne-1
      medianCI <- 2*medianCIZeroOne-1
      upperCI <- 2*upperCIZeroOne-1
    }
  }
  result <- c(lowerCI, medianCI, upperCI)
  return(result)
}

