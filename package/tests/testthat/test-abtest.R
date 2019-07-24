context('BFDA for abtest')

# Check sampling function
test_that("AB sampling works", {
  
  # Draw sample for Odds Ratio and log Odds Ratio
  set.seed(1234)
  SAMP1  <- BFDA:::sample.abtest(10000, 2, options.sample = list(effecttype = "OR"))
  empES1 <- BFDA:::freq.test.abtest(SAMP1, alternative = "two.sided", options.sample = list(effecttype = "OR"))$emp.ES
  set.seed(1234)
  SAMP2  <- BFDA:::sample.abtest(10000, log(2), options.sample = list(effecttype = "logOR"))
  empES2 <- BFDA:::freq.test.abtest(SAMP2, alternative = "two.sided", options.sample = list(effecttype = "logOR"))$emp.ES
  
  expect_identical(SAMP1, SAMP2) # sampling log odds ratio and log(odds ratio) should be same 
  expect_equal(empES1, 2, tolerance = 0.05)
  expect_equal(empES2, log(2), tolerance = 0.05)
  
  # Draw sample for Relative Risk
  set.seed(1234)
  SAMP3  <- BFDA:::sample.abtest(10000, 2, options.sample = list(effecttype = "RR"))
  empES3 <- BFDA:::freq.test.abtest(SAMP3, alternative = "two.sided", options.sample = list(effecttype = "RR"))$emp.ES
  
  expect_equal(empES3, 2, tolerance = 0.05)
  
  # Draw sample for Absolute Risk
  set.seed(1234)
  SAMP4 <- BFDA:::sample.abtest(10000, 0.5, options.sample = list(effecttype = "AR"))
  empES3 <- BFDA:::freq.test.abtest(SAMP4, alternative = "two.sided", options.sample = list(effecttype = "AR"))$emp.ES
  
  expect_equal(empES3, 0.5, tolerance = 0.05)
  expect_error(BFDA:::sample.abtest(10000, 2, options.sample = list(effecttype = "AR")),
               "Absolute risk can only take values between 0 and 1")
  
})
# Check internal functions for Bayesian AB test

test_that("AB test works", {
  
  set.seed(1234)
  SAMP <- BFDA:::sample.abtest(100, 2, options.sample = list(effecttype = "OR"))
  freq.test <- BFDA:::freq.test.abtest(SAMP, alternative="two.sided", options.sample = list(effecttype = "OR"))
  
  y1 <- colSums(SAMP)[1]
  y2 <- colSums(SAMP)[2]
  n1 <- n2 <- nrow(SAMP)
  data <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2)
  prior_par <- list(mu_psi = 0, sigma_psi = 1, mu_beta = 0, sigma_beta = 1)
  prior <- list("normal", list(prior.mean = 0, prior.variance = 1))
  
  # Alternative = two.sided
  # Compare internal functions
  
  bf1.1 <- abtest:::ab_test(data = data, prior_par = prior_par)$bf
  bf1.2 <- exp(BFDA:::BF.test.abtest(SAMP = SAMP, alternative="two.sided", freq.test=freq.test, prior=list("normal", list(prior.mean = 0, prior.variance = 1))))
  
  expect_equal(bf1.1$bf10, bf1.2, tolerance = 0.01) # internal functions doing the same?
  
  # Alternative = greater
  # Compare internal functions
  
  bf2.1 <- bf1.1$bfplus0
  bf2.2 <- exp(BFDA:::BF.test.abtest(SAMP = SAMP, alternative="greater", freq.test=freq.test, prior=list("normal", list(prior.mean = 0, prior.variance = 1))))
  
  expect_equal(bf2.1, bf2.2, tolerance = 0.01) # internal functions doing the same?
  
  # Alternative = less
  # Compare internal functions
  
  bf3.1 <- bf1.1$bfminus0
  bf3.2 <- exp(BFDA:::BF.test.abtest(SAMP = SAMP, alternative="less", freq.test=freq.test, prior=list("normal", list(prior.mean = 0, prior.variance = 1))))
  
  expect_equal(bf3.1, bf3.2, tolerance = 0.01) # internal functions doing the same?
  
})

test_that("BFDA for AB test works", {
  
  # fixed.n BFDA
  
  BFDA.res1.1 <- BFDA.sim(expected.ES=2, type="abtest", prior = list("normal", list(prior.mean = 0, prior.variance = 1)), n.max=100, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative="two.sided", verbose=FALSE, cores=1, ETA=FALSE, options.sample=list(effecttype = "OR"), seed=1234)
  BFDA.res1.2 <- BFDA.sim(expected.ES=2, type="abtest", n.max=100, design="fixed.n", boundary=Inf, B=2, stepsize=NA, alternative="two.sided", verbose=FALSE, cores=1, ETA=FALSE, options.sample=list(effecttype = "OR"), seed=1234)
  
  expect_equal(BFDA.res1.1$sim, BFDA.res1.2$sim, tolerance = 0.01) # because default is applied
  
  expect_warning(BFDA.sim(expected.ES=2, type="abtest", prior = list("normal", list(prior.location = 0, prior.variance = 1)), n.max=100, design=c("fixed.n"), boundary=Inf, B=2, stepsize=NA, alternative="two.sided", verbose=FALSE, cores=1, ETA=FALSE, options.sample=list(effecttype = "OR"), seed=1234),
                 "Normal distribution only takes parameters prior.mean and prior.variance.")
  expect_warning(BFDA.sim(expected.ES=2, type="abtest", prior = list("normal", list(prior.location = 0, prior.variance = 1)), n.max=100, design=c("fixed.n"), boundary=Inf, B=2, stepsize=NA, alternative="two.sided", verbose=FALSE, cores=1, ETA=FALSE, options.sample=list(effecttype = "OR"), seed=1234),
                 "Prior mean not defined. Default specification will be used.")
  expect_warning(BFDA.sim(expected.ES=2, type="abtest", prior = list("normal", list(prior.mean = 0)), n.max=100, design=c("fixed.n"), boundary=Inf, B=2, stepsize=NA, alternative="two.sided", verbose=FALSE, cores=1, ETA=FALSE, options.sample=list(effecttype = "OR"), seed=1234),
                 "Prior variance not defined. Default specification will be used.")
  
  # Sequential BFDA
  
  BFDA.res2.1 <- BFDA.sim(expected.ES=2, type="abtest", prior = list("normal", list(prior.mean = 0, prior.variance = 1)), n.min = 10, n.max=20, design="sequential", boundary=6, B=2, stepsize=NA, alternative="two.sided", verbose=FALSE, cores=1, ETA=FALSE, options.sample=list(effecttype = "OR"), seed=1234)
  BFDA.res2.2 <- BFDA.sim(expected.ES=2, type="abtest", n.min = 10, n.max=20, design="sequential", boundary=6, B=2, stepsize=NA, alternative="two.sided", verbose=FALSE, cores=1, ETA=FALSE, options.sample=list(effecttype = "OR"), seed=1234)
  
  expect_equal(BFDA.res2.1$sim, BFDA.res2.2$sim, tolerance = 0.01) # because default is applied
  
  expect_error(BFDA.sim(expected.ES=2, type="abtest", prior = 1, n.min = 10, n.max=20, design="sequential", boundary=6, B=2, stepsize=NA, alternative="two.sided", verbose=FALSE, cores=1, ETA=FALSE, options.sample=list(effecttype = "OR"), seed=1234),
               "Argument prior needs to be specified as a list.")
    
})
