context('BFDA for correlation')

# Check internal functions for Bayesian correlation test

test_that("Correlation test works", {
  
  set.seed(1234)
  SAMP <- BFDA:::getBiCop(20, 0.3)
  freq.test <- BFDA:::freq.test.correlation(SAMP, alternative="two.sided")
  
  # Alternative = two.sided
  # Compares results of BFDA functions and BayesFactor package 

  bf1.1 <- exp(BayesFactor::correlationBF(SAMP[,1], SAMP[,2], nullInterval=NULL, posterior=FALSE, rscale = 1)@bayesFactor$bf)
  bf1.2 <- BFDA:::bfPearsonCorrelation(nrow(SAMP), as.numeric(freq.test$emp.ES), kappa=1)$bf10
  bf1.3 <- exp(BFDA:::BF.test.correlation(SAMP = SAMP, alternative="two.sided", freq.test=freq.test, prior=list("stretchedbeta", list(prior.kappa = 1))))
  
  expect_equal(bf1.1, bf1.2, tolerance = 0.01) # BayesFactor vs. BFDA
  expect_equal(bf1.2, bf1.3, tolerance = 0.01) # BFDA internal functions
  
  # Alternative = greater
  # Compares results of BFDA functions and BayesFactor package
  
  bf2.1 <- exp(BayesFactor::correlationBF(SAMP[,1], SAMP[,2], nullInterval=c(0, 1), posterior=FALSE, rscale = 1)@bayesFactor$bf[1])
  bf2.2 <- BFDA:::bfPearsonCorrelation(nrow(SAMP), as.numeric(freq.test$emp.ES), kappa=1)$bfPlus0
  bf2.3 <- exp(BFDA:::BF.test.correlation(SAMP = SAMP, alternative="greater", freq.test=freq.test, prior=list("stretchedbeta", list(prior.kappa = 1))))
  
  expect_equal(bf2.1, bf2.2, tolerance = 0.01) # BayesFactor vs. BFDA
  expect_equal(bf2.2, bf2.3, tolerance = 0.01) # BFDA internal functions
  
  # Alternative = less
  # Compares results of BFDA functions and BayesFactor package
  
  bf3.1 <- exp(BayesFactor::correlationBF(SAMP[,1], SAMP[,2], nullInterval=c(-1, 0), posterior=FALSE, rscale = 1)@bayesFactor$bf[1])
  bf3.2 <- BFDA:::bfPearsonCorrelation(nrow(SAMP), as.numeric(freq.test$emp.ES), kappa=1)$bfMin0
  bf3.3 <- exp(BFDA:::BF.test.correlation(SAMP = SAMP, alternative="less", freq.test=freq.test, prior=list("stretchedbeta", list(prior.kappa = 1))))
  
  expect_equal(bf3.1, bf3.2, tolerance = 0.01) # BayesFactor vs. BFDA
  expect_equal(bf3.2, bf3.3, tolerance = 0.01) # BFDA internal functions
  
})

# Check the BFDA for the correlation test

test_that("BFDA for correlation test works", {
  
  # fixed.n BFDA
  
  BFDA.res1.1 <- BFDA.sim(expected.ES=0.3, type="correlation", prior = list("stretchedbeta", list(prior.kappa = 1)), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  BFDA.res1.2 <- BFDA.sim(expected.ES=0.3, type="correlation", n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1) #should apply default
  
  expect_equal(BFDA.res1.1, BFDA.res1.2, tolerance = 0.01) # because default is applied
  
  expect_warning(BFDA.sim(expected.ES=0.3, type="correlation", prior = list("stretchedbeta", list(prior.beta = 1)), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Stretched beta prior only takes parameter prior.kappa. Default value will be applied.")
  
  expect_error(BFDA.sim(expected.ES=0.3, type="correlation", n.max=20, prior = 1, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("less"), verbose=FALSE, cores=1),
               "Argument prior needs to be specified as a list.")
  
  # sequential BFDA
  
  BFDA.res2.1 <- BFDA.sim(expected.ES=0.3, type="correlation", prior = list("stretchedbeta", list(prior.kappa = 1)), n.max=20, design="sequential", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  BFDA.res2.2 <- BFDA.sim(expected.ES=0.3, type="correlation", n.max=20, design="sequential", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1) #should apply default
  
  expect_equal(BFDA.res2.1, BFDA.res2.2, tolerance = 0.01) # because default is applied
  
  expect_warning(BFDA.sim(expected.ES=0.3, type="correlation", prior = list("stretchedbeta", list(prior.beta = 1)), n.max=20, design="sequential", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Stretched beta prior only takes parameter prior.kappa. Default value will be applied.")
  
  expect_error(BFDA.sim(expected.ES=0.3, type="correlation", n.max=20, prior = 1, design="sequential", boundary=Inf, B=2, stepsize=NA, alternative=c("less"), verbose=FALSE, cores=1),
               "Argument prior needs to be specified as a list.")
  
})

