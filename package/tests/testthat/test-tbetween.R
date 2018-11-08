context('BFDA for t.between')

# Check internal functions for Bayesian t-test

test_that("Independent samples t-test works", {
  
  set.seed(1234)
  SAMP <- BFDA:::sample.t.between(20, 0.3)
  freq.test <- BFDA:::freq.test.t.between(SAMP, alternative="two.sided")
  
  # Alternative = two.sided, default Cauchy Prior
  # Compares results of BFDA functions and BayesFactor package 
  
  bf1.1 <- as.numeric(BayesFactor::ttest.tstat(freq.test$statistic, n1=nrow(SAMP), n2=nrow(SAMP), nullInterval=NULL, simple=TRUE))
  bf1.2 <- BFDA:::bf10_t(t = as.numeric(freq.test$statistic), n1 = nrow(SAMP), n2=nrow(SAMP), independentSamples = TRUE, prior.location = 0, prior.scale = sqrt(2)/2, prior.df = 1)$BF10
  bf1.3 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="two.sided", freq.test=freq.test, prior=list("Cauchy", list(prior.location = 0, prior.scale = sqrt(2)/2))))
  
  expect_equal(bf1.1, bf1.2, tolerance = 0.01)
  expect_equal(bf1.2, bf1.3, tolerance = 0.01)
  
  # Alternative = greater, default Cauchy Prior
  # Compares results of BFDA functions and BayesFactor package 
  
  bf2.1 <- as.numeric(BayesFactor::ttest.tstat(freq.test$statistic, n1=nrow(SAMP), n2=nrow(SAMP), nullInterval=c(0, Inf), simple=TRUE))
  bf2.2 <- BFDA:::bf10_t(t = as.numeric(freq.test$statistic), n1 = nrow(SAMP), n2=nrow(SAMP), independentSamples = TRUE, prior.location = 0, prior.scale = sqrt(2)/2, prior.df = 1)$BFplus0
  bf2.3 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="greater", freq.test=freq.test, prior=list("Cauchy", list(prior.location = 0, prior.scale = sqrt(2)/2))))
  
  expect_equal(bf2.1, bf2.2, tolerance = 0.01)
  expect_equal(bf2.2, bf2.3, tolerance = 0.01)
  
  # Alternative = less, default Cauchy Prior
  # Compares results of BFDA functions and BayesFactor package 
  
  bf3.1 <- as.numeric(BayesFactor::ttest.tstat(freq.test$statistic, n1=nrow(SAMP), n2=nrow(SAMP), nullInterval=c(-Inf, 0), simple=TRUE))
  bf3.2 <- BFDA:::bf10_t(t = as.numeric(freq.test$statistic), n1 = nrow(SAMP), n2=nrow(SAMP), independentSamples = TRUE, prior.location = 0, prior.scale = sqrt(2)/2, prior.df = 1)$BFmin0
  bf3.3 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="less", freq.test=freq.test, prior=list("Cauchy", list(prior.location = 0, prior.scale = sqrt(2)/2))))
  
  expect_equal(bf3.1, bf3.2, tolerance = 0.01)
  expect_equal(bf3.2, bf3.3, tolerance = 0.01)
  
  # t Prior
  # Compares results of t priors to Cauchy priors (t with 1 df = Cauchy)
  
  bf4.1 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="two.sided", freq.test=freq.test, prior=list("t", list(prior.location = 0, prior.scale = sqrt(2)/2, prior.df = 1))))
  bf4.2 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="greater", freq.test=freq.test, prior=list("t", list(prior.location = 0, prior.scale = sqrt(2)/2, prior.df = 1))))
  bf4.3 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="less", freq.test=freq.test, prior=list("t", list(prior.location = 0, prior.scale = sqrt(2)/2, prior.df = 1))))
  
  expect_equal(bf1.2, bf4.1, tolerance = 0.01)
  expect_equal(bf2.2, bf4.2, tolerance = 0.01)
  expect_equal(bf3.2, bf4.3, tolerance = 0.01)

  # Normal Prior
  # Compares results of normal priors to t priors (t distribution with many df -> normal distribution)
  
  bf5.1.1 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="two.sided", freq.test=freq.test, prior=list("t", list(prior.location = 0, prior.scale = 1, prior.df = 10000))))
  bf5.2.1 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="greater", freq.test=freq.test, prior=list("t", list(prior.location = 0, prior.scale = 1, prior.df = 10000))))
  bf5.3.1 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="less", freq.test=freq.test, prior=list("t", list(prior.location = 0, prior.scale = 1, prior.df = 100000))))
  
  bf5.1.2 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="two.sided", freq.test=freq.test, prior=list("normal", list(prior.mean = 0, prior.variance = 1))))
  bf5.2.2 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="greater", freq.test=freq.test, prior=list("normal", list(prior.mean = 0, prior.variance = 1))))
  bf5.3.2 <- exp(BFDA:::BF.test.t.between(SAMP = SAMP, alternative="less", freq.test=freq.test, prior=list("normal", list(prior.mean = 0, prior.variance = 1))))
  
  expect_equal(bf5.1.1, bf5.1.2, tolerance = 0.01)
  expect_equal(bf5.2.1, bf5.2.2, tolerance = 0.01)
  expect_equal(bf5.3.1, bf5.3.2, tolerance = 0.01)
  
})

# Check the BFDA for the independent samples t-test

test_that("BFDA for independent samples t-test works", {
  
  # fixed.n BFDA with Cauchy prior
  
  BFDA.res1.1 <- BFDA.sim(expected.ES=0.3, type="t.between", prior = list("Cauchy", list(prior.location = 0, prior.scale = sqrt(2)/2)), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  BFDA.res1.2 <- BFDA.sim(expected.ES=0.3, type="t.between", n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1) #should apply defaults

  expect_equal(BFDA.res1.1$sim, BFDA.res1.2$sim, tolerance = 0.01) #because defaults are applied in 1.2
  
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("Cauchy", list(prior.location = 0)), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Prior scale not defined. Default specification will be used.")
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("Cauchy", list()), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Prior location not defined. Default specification will be used.")
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("Cauchy", list()), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Prior scale not defined. Default specification will be used.")
  
  expect_error(BFDA.sim(expected.ES=0.3, type="t.between", prior = "something", n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
               "Argument prior needs to be specified as a list.")
  
  # fixed.n BFDA with t prior
  
  BFDA.res2.1 <- BFDA.sim(expected.ES=0.3, type="t.between", prior = list("t", list(prior.location = 0, prior.scale = sqrt(2)/2, prior.df = 1)), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  BFDA.res2.2 <- BFDA.sim(expected.ES=0.3, type="t.between", n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  
  
  expect_equal(BFDA.res2.1$sim, BFDA.res1.1$sim, tolerance = 0.01) # because t with 1 df = Cauchy
  expect_equal(BFDA.res2.1$sim, BFDA.res2.2$sim, tolerance = 0.01) # because defaults are applied
  
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("t", list()), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Prior location not defined. Default specification will be used.")
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("t", list()), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Prior scale not defined. Default specification will be used.")
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("t", list()), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Prior degrees of freedom not defined. Default specifications will be used.")
  
  # fixed.n BFDA with normal prior
  
  BFDA.res3.1 <- BFDA.sim(expected.ES=0.3, type="t.between", prior = list("t", list(prior.location = 0, prior.scale = 1, prior.df = 100000)), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  BFDA.res3.2 <- BFDA.sim(expected.ES=0.3, type="t.between", prior = list("normal", list(prior.mean = 0, prior.variance = 1)), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  BFDA.res3.3 <- suppressWarnings(BFDA.sim(expected.ES=0.3, type="t.between", n.max=20, list("normal", list()), design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1))
  
  
  expect_equal(BFDA.res3.1$sim, BFDA.res3.2$sim, tolerance = 0.01) # because t with many df -> normal
  expect_equal(BFDA.res3.2$sim, BFDA.res3.3$sim, tolerance = 0.01) # because defaults are applied
  
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("normal", list(prior.mean = 0)), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Prior variance not defined. Default specification will be used.")
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("normal", list()), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Prior mean not defined. Default specification will be used.")
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("normal", list()), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Prior variance not defined. Default specification will be used.")
  expect_warning(BFDA.sim(expected.ES=0.3, type="t.between", prior = list("normal", list(prior.location = 0)), n.max=20, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1),
                 "Normal distribution only takes parameters prior.mean and prior.variance.")
  
  # sequential BFDA with Cauchy prior
  
  BFDA.res4.1 <- BFDA.sim(expected.ES=0.3, type="t.between", prior = list("Cauchy", list(prior.location = 0, prior.scale = sqrt(2)/2)), n.min = 10, n.max=20, design="sequential", boundary=6, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  BFDA.res4.2 <- BFDA.sim(expected.ES=0.3, type="t.between", n.min = 10, n.max=20, design="sequential", boundary=6, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1) #should apply defaults
  
  expect_equal(BFDA.res4.1$sim, BFDA.res4.2$sim, tolerance = 0.01) #because defaults are applied in 1.2
  
  # sequential BFDA with t prior
  
  BFDA.res5.1 <- BFDA.sim(expected.ES=0.3, type="t.between", prior = list("t", list(prior.location = 0, prior.scale = sqrt(2)/2, prior.df = 1)), n.min = 10, n.max=20, design="sequential", boundary=6, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  BFDA.res5.2 <- BFDA.sim(expected.ES=0.3, type="t.between", n.min = 10, n.max=20, design="sequential", boundary=6, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1) #should apply defaults
  
  expect_equal(BFDA.res5.1$sim, BFDA.res5.2$sim, tolerance = 0.01) #because defaults are applied in 1.2
  
  # sequential BFDA with normal prior
  
  BFDA.res6.1 <- BFDA.sim(expected.ES=0.3, type="t.between", prior = list("normal", list(prior.mean = 0, prior.variance = 1)), n.min = 10, n.max=20, design="sequential", boundary=6, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  BFDA.res6.2 <- BFDA.sim(expected.ES=0.3, type="t.between", prior = list("t", list(prior.location = 0, prior.scale = 1, prior.df = 100000)), n.min = 10, n.max=20, design="sequential", boundary=6, B=2, stepsize=NA, alternative=c("two.sided"), verbose=FALSE, cores=1)
  
  expect_equal(BFDA.res6.1$sim, BFDA.res6.2$sim, tolerance = 0.01) #because defaults are applied in 1.2
  
})

