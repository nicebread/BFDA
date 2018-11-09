## =============================================================================
## These are the the functions for the AB-test
## =============================================================================
#  The functions are a simplified version of the abtest R package by 
#  Quentin Gronau

# ------------------------------------------------------------------------------
# Helper functions

# gradient of minus l0(beta) ----
gradient_minl0 <- function(par, data, mu_beta = 0, sigma_beta = 1) {
  
  # parameter
  beta <- par[1]
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  out <- - (y1 + y2 - (n1 + n2 - y1 - y2) * exp(beta)) /
    (1 + exp(beta)) + (beta - mu_beta) / sigma_beta^2
  
  return(out)
  
}

# minus l0(beta) ----
minl0 <- function(par, data, mu_beta = 0, sigma_beta = 1) {
  
  # parameter
  beta <- par[1]
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  # log probabilities
  logp <- beta - VGAM::log1pexp(beta)
  log1minp <- - VGAM::log1pexp(beta)
  
  out <- - (y1 + y2) * logp - (n1 + n2 - y1 - y2) * log1minp -
    stats::dnorm(beta, mu_beta, sigma_beta, log = TRUE)
  
  return(out)
  
}

# gradient of minus l0(beta) -----
gradient_minl0 <- function(par, data, mu_beta = 0,
                           sigma_beta = 1) {
  
  # parameter
  beta <- par[1]
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  out <- - (y1 + y2 - (n1 + n2 - y1 - y2) * exp(beta)) /
    (1 + exp(beta)) + (beta - mu_beta) / sigma_beta^2
  
  return(out)
  
}

# gradient of minus l(beta, psi) ----

gradient_minl <- function(par, data, mu_beta = 0, sigma_beta = 1, mu_psi = 0, sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  psi <- par[2]
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  partial_beta <- (y1 - (n1 - y1) * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi)) + (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi)) - (beta - mu_beta) / sigma_beta^2
  
  partial_psi <- .5 * (((n1 - y1) * exp(beta - .5 * psi) - y1) /
                         (1 + exp(beta - .5 * psi)) +
                         (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
                         (1 + exp(beta + .5 * psi))) -
    (psi - mu_psi) / sigma_psi^2
  
  return(- c(partial_beta, partial_psi))
  
}

# minus l(beta, psi) -----

minl <- function(par, data, mu_beta = 0, sigma_beta = 1,
                 mu_psi = 0, sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  psi <- par[2]
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  # log probabilities
  logp1 <- beta - .5 * psi - VGAM::log1pexp(beta - .5 * psi)
  log1minp1 <- - VGAM::log1pexp(beta - .5 * psi)
  logp2 <- beta + .5 * psi - VGAM::log1pexp(beta + .5 * psi)
  log1minp2 <- - VGAM::log1pexp(beta + .5 * psi)
  
  out <- y1 * logp1 + (n1 - y1) * log1minp1 +
    y2 * logp2 + (n2 - y2) * log1minp2 +
    stats::dnorm(beta, mu_beta, sigma_beta, log = TRUE) +
    stats::dnorm(psi, mu_psi, sigma_psi, log = TRUE)
  
  return(- out)
  
}

# inverse of Hessian for minus l0(beta) -----
inverse_hessian_minl0 <- function(par, data, sigma_beta = 1) {
  
  hessian <- hessian_minl0(par = par, data = data, sigma_beta = sigma_beta)
  return(1 / hessian)
  
}

# Hessian for minus l0(beta) -----
hessian_minl0 <- function(par, data, sigma_beta = 1) {
  
  # parameter
  beta <- par[1]
  
  # data
  n1 <- data$n1
  n2 <- data$n2
  
  hessian <- (n1 + n2) * exp(beta) / (1 + exp(beta))^2 + 1 / sigma_beta^2
  return(hessian)
  
}

# inverse of Hessian for minus l(beta, psi) ----
inverse_hessian_minl <- function(par, data, sigma_beta = 1,
                                 sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  psi <- par[2]
  
  # data
  n1 <- data$n1
  n2 <- data$n2
  
  det_hessian <- det_hessian_minl(par = par, data = data,
                                  sigma_beta = sigma_beta,
                                  sigma_psi = sigma_beta)
  
  partial_beta2 <-  (n1 * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi))^2 + (n2 * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi))^2 + 1 / sigma_beta^2
  
  partial_psi2 <- .25 * ((n1 * exp(beta - .5 * psi)) /
                           (1 + exp(beta - .5 * psi))^2 +
                           (n2 * exp(beta + .5 * psi)) /
                           (1 + exp(beta + .5 * psi))^2) +
    1 / sigma_psi^2
  
  minus_partial_beta_psi <- .5 * ((n1 * exp(beta - .5 * psi)) /
                                    (1 + exp(beta - .5 * psi))^2 -
                                    (n2 * exp(beta + .5 * psi)) /
                                    (1 + exp(beta + .5 * psi))^2)
  
  inv_hessian <- 1 / det_hessian *
    matrix(c(partial_psi2, rep(minus_partial_beta_psi, 2),
             partial_beta2), 2, 2, byrow = TRUE)
  
  return(inv_hessian)
  
}

# determinant of Hessian for minus l(beta, psi) -----
det_hessian_minl <- function(par, data, sigma_beta = 1, sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  psi <- par[2]
  
  # data
  n1 <- data$n1
  n2 <- data$n2
  
  hessian <- hessian_minl(par = par, data = data,
                          sigma_beta = sigma_beta,
                          sigma_psi = sigma_psi)
  det_hessian <- hessian[1,1] * hessian[2,2] - hessian[1,2]^2
  
  return(det_hessian)
  
}

# Hessian for minus l(beta, psi) -----
hessian_minl <- function(par, data, sigma_beta = 1, sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  psi <- par[2]
  
  # data
  n1 <- data$n1
  n2 <- data$n2
  
  partial_beta2 <-  (n1 * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi))^2 + (n2 * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi))^2 + 1 / sigma_beta^2
  
  partial_psi2 <- .25 * ((n1 * exp(beta - .5 * psi)) /
                           (1 + exp(beta - .5 * psi))^2 +
                           (n2 * exp(beta + .5 * psi)) /
                           (1 + exp(beta + .5 * psi))^2) +
    1 / sigma_psi^2
  
  partial_beta_psi <- - .5 * ((n1 * exp(beta - .5 * psi)) /
                                (1 + exp(beta - .5 * psi))^2 -
                                (n2 * exp(beta + .5 * psi)) /
                                (1 + exp(beta + .5 * psi))^2)
  
  hessian <- matrix(c(partial_beta2, rep(partial_beta_psi, 2),
                      partial_psi2), 2, 2, byrow = TRUE)
  return(hessian)
  
}

# minus l+(beta, eta)
minlplus <- function(par, data, mu_beta = 0, sigma_beta = 1,
                     mu_psi = 0, sigma_psi = 1) {
  
  
  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- exp(xi)
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  # log probabilities
  logp1 <- beta - .5 * psi - log1pexp(beta - .5 * psi)
  log1minp1 <- - log1pexp(beta - .5 * psi)
  logp2 <- beta + .5 * psi - log1pexp(beta + .5 * psi)
  log1minp2 <- - log1pexp(beta + .5 * psi)
  
  out <- y1 * logp1 + (n1 - y1) * log1minp1 +
    y2 * logp2 + (n2 - y2) * log1minp2 +
    dnorm(beta, mu_beta, sigma_beta, log = TRUE) +
    dnorm(psi, mu_psi, sigma_psi, log = TRUE) -
    pnorm(0, mu_psi, sigma_psi, lower.tail = FALSE, log.p = TRUE) + xi
  
  return(- out)
  
}

# gradient of minus l+(beta, xi)
gradient_minlplus <- function(par, data, mu_beta = 0,
                              sigma_beta = 1, mu_psi = 0,
                              sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- exp(xi)
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  partial_beta <- (y1 - (n1 - y1) * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi)) + (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi)) - (beta - mu_beta) / sigma_beta^2
  
  partial_xi <- .5 * psi * (((n1 - y1) * exp(beta - .5 * psi) - y1) /
                              (1 + exp(beta - .5 * psi)) +
                              (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
                              (1 + exp(beta + .5 * psi))) -
    psi * (psi - mu_psi) / sigma_psi^2 + 1
  
  return(- c(partial_beta, partial_xi))
  
}

# Hessian for minus l+(beta, xi)
hessian_minlplus <- function(par, data, sigma_beta = 1,
                             mu_psi = 0, sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- exp(xi)
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  partial_beta2 <-  (n1 * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi))^2 + (n2 * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi))^2 + 1 / sigma_beta^2
  
  partial_xi2 <- - .5 * psi * (((n1 - y1) * exp(beta - .5 * psi) - y1) /
                                 (1 + exp(beta - .5 * psi)) +
                                 (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
                                 (1 + exp(beta + .5 * psi)) -
                                 .5 * psi * n1 * exp(beta - .5 * psi) /
                                 (1 + exp(beta - .5 * psi))^2 -
                                 .5 * psi * n2 * exp(beta + .5 * psi) /
                                 (1 + exp(beta + .5 * psi))^2) +
    psi * (2 * psi - mu_psi) / sigma_psi^2
  
  partial_beta_xi <- - .5 * psi * ((n1 * exp(beta - .5 * psi)) /
                                     (1 + exp(beta - .5 * psi))^2 -
                                     (n2 * exp(beta + .5 * psi)) /
                                     (1 + exp(beta + .5 * psi))^2)
  
  hessian <- matrix(c(partial_beta2, rep(partial_beta_xi, 2),
                      partial_xi2), 2, 2, byrow = TRUE)
  
  return(hessian)
  
}

# inverse of Hessian for minus l+(beta, xi)
inverse_hessian_minlplus <- function(par, data, sigma_beta = 1, mu_psi = 0,
                                     sigma_psi = 1) {
  
  
  det_hessian <- det_hessian_minlplus(par = par, data = data,
                                      sigma_beta = sigma_beta,
                                      mu_psi = mu_psi,
                                      sigma_psi = sigma_beta)
  
  hessian <- hessian_minlplus(par = par, data = data,
                              sigma_beta = sigma_beta,
                              mu_psi = mu_psi,
                              sigma_psi = sigma_beta)
  
  inv_hessian <- 1 / det_hessian *
    matrix(c(hessian[2,2], rep(- hessian[1,2], 2),
             hessian[1,1]), 2, 2, byrow = TRUE)
  
  return(inv_hessian)
  
}

# minus l-(beta, xi)
minlminus <- function(par, data, mu_beta = 0, sigma_beta = 1,
                      mu_psi = 0, sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- - exp(xi)
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  # log probabilities
  logp1 <- beta - .5 * psi - log1pexp(beta - .5 * psi)
  log1minp1 <- - log1pexp(beta - .5 * psi)
  logp2 <- beta + .5 * psi - log1pexp(beta + .5 * psi)
  log1minp2 <- - log1pexp(beta + .5 * psi)
  
  out <- y1 * logp1 + (n1 - y1) * log1minp1 +
    y2 * logp2 + (n2 - y2) * log1minp2 +
    dnorm(beta, mu_beta, sigma_beta, log = TRUE) +
    dnorm(psi, mu_psi, sigma_psi, log = TRUE) -
    pnorm(0, mu_psi, sigma_psi, log.p = TRUE) + xi
  
  return(- out)
  
}

# gradient of minus l+(beta, xi)
gradient_minlminus <- function(par, data, mu_beta = 0,
                               sigma_beta = 1, mu_psi = 0,
                               sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- - exp(xi)
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  partial_beta <- (y1 - (n1 - y1) * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi)) + (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi)) - (beta - mu_beta) / sigma_beta^2
  
  partial_xi <- - .5 * psi * ((y1 - (n1 - y1) * exp(beta - .5 * psi)) /
                                (1 + exp(beta - .5 * psi)) +
                                ((n2 - y2) * exp(beta + .5 * psi) - y2) /
                                (1 + exp(beta + .5 * psi))) -
    psi * (psi - mu_psi) / sigma_psi^2 + 1
  
  return(- c(partial_beta, partial_xi))
  
}

# inverse of Hessian for minus l-(beta, xi)
inverse_hessian_minlminus <- function(par, data, sigma_beta = 1,
                                      mu_psi = 0, sigma_psi = 1) {
  
  det_hessian <- det_hessian_minlminus(par = par, data = data,
                                       sigma_beta = sigma_beta,
                                       mu_psi = mu_psi,
                                       sigma_psi = sigma_beta)
  
  hessian <- hessian_minlminus(par = par, data = data,
                               sigma_beta = sigma_beta,
                               mu_psi = mu_psi,
                               sigma_psi = sigma_beta)
  
  inv_hessian <- 1 / det_hessian *
    matrix(c(hessian[2,2], rep(- hessian[1,2], 2),
             hessian[1,1]), 2, 2, byrow = TRUE)
  
  return(inv_hessian)
  
}

# determinant of Hessian for minus l-(beta, xi)
det_hessian_minlminus <- function(par, data, sigma_beta = 1, mu_psi = 0,
                                  sigma_psi = 1) {
  
  hessian <- hessian_minlminus(par = par, data = data,
                               sigma_beta = sigma_beta,
                               mu_psi = mu_psi,
                               sigma_psi = sigma_psi)
  
  det_hessian <- hessian[1,1] * hessian[2,2] - hessian[1,2]^2
  
  return(det_hessian)
  
}

# Hessian for minus l-(beta, xi)
hessian_minlminus <- function(par, data, sigma_beta = 1,
                              mu_psi = 0, sigma_psi = 1) {
  
  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- - exp(xi)
  
  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2
  
  partial_beta2 <-  (n1 * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi))^2 + (n2 * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi))^2 + 1 / sigma_beta^2
  
  partial_xi2 <- - .5 * psi * (((n1 - y1) * exp(beta - .5 * psi) - y1) /
                                 (1 + exp(beta - .5 * psi)) +
                                 (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
                                 (1 + exp(beta + .5 * psi)) -
                                 .5 * psi * n1 * exp(beta - .5 * psi) /
                                 (1 + exp(beta - .5 * psi))^2 -
                                 .5 * psi * n2 * exp(beta + .5 * psi) /
                                 (1 + exp(beta + .5 * psi))^2) +
    psi * (2 * psi - mu_psi) / sigma_psi^2
  
  partial_beta_xi <- - .5 * psi * ((n1 * exp(beta - .5 * psi)) /
                                     (1 + exp(beta - .5 * psi))^2 -
                                     (n2 * exp(beta + .5 * psi)) /
                                     (1 + exp(beta + .5 * psi))^2)
  
  hessian <- matrix(c(partial_beta2, rep(partial_beta_xi, 2),
                      partial_xi2), 2, 2, byrow = TRUE)
  
  return(hessian)
  
}

# determinant of Hessian for minus l+(beta, xi)
det_hessian_minlplus <- function(par, data, sigma_beta = 1, mu_psi = 0,
                                 sigma_psi = 1) {
  
  hessian <- hessian_minlplus(par = par, data = data,
                              sigma_beta = sigma_beta,
                              mu_psi = mu_psi,
                              sigma_psi = sigma_psi)
  det_hessian <- hessian[1,1] * hessian[2,2] - hessian[1,2]^2
  
  return(det_hessian)
  
}

# ------------------------------------------------------------------------------
# Function for conducting a Bayesian A/B test.

#' @importFrom Rcpp sourceCpp
#' @importFrom stats dnorm median optim pnorm qlogis
#' @importFrom VGAM log1pexp
#' @importFrom mvtnorm dmvt rmvt
#' @importFrom sn psn selm
#' @useDynLib BFDA 

ab_test <- function(data, prior_par, nsamples = 1e4, is_df = 5) {
  
  # extract prior parameters
  mu_psi <- prior_par[["mu_psi"]]
  sigma_psi <- prior_par[["sigma_psi"]]
  mu_beta <- prior_par[["mu_beta"]]
  sigma_beta <- prior_par[["sigma_beta"]]
  
  ### H1 and H0 ###
  
  # starting values
  start_beta <- stats::qlogis((data$y1 + data$y2) / (data$n1 + data$n2))
  eta1 <- stats::qlogis(data$y1 / data$n1)
  eta2 <- stats::qlogis(data$y2 / data$n2)
  
  if (is.infinite(start_beta)) {
    start_beta <- 0 # a finite number
  }
  if (is.infinite(eta1)) {
    eta1 <- 0 # a finite number
  }
  if (is.infinite(eta2)) {
    eta2 <- 0 # a finite number
  }
  
  start_psi <- eta2 - eta1
  
  # optimize
  o0 <- stats::optim(par = start_beta, fn = minl0, gr = gradient_minl0,
              data = data, method = "BFGS", mu_beta = mu_beta,
              sigma_beta = sigma_beta)
  o1 <- stats::optim(par = c(start_beta, start_psi), fn = minl,
              gr = gradient_minl, data = data, method = "BFGS",
              mu_beta = mu_beta, sigma_beta = sigma_beta,
              mu_psi = mu_psi, sigma_psi = sigma_psi)
  
  # approximate posterior modes
  mode_H0 <- o0$par
  names(mode_H0) <- "beta"
  mode_H1 <- o1$par
  parnames <-  c("beta", "psi")
  names(mode_H1) <- parnames
  
  # approximate posterior variance/covariance
  var_H0 <- inverse_hessian_minl0(par = o0$par, data = data,
                                  sigma_beta = sigma_beta)
  names(var_H0) <- "beta"
  cov_H1 <- inverse_hessian_minl(par = o1$par, data = data,
                                 sigma_beta = sigma_beta,
                                 sigma_psi = sigma_psi)
  colnames(cov_H1) <- parnames
  rownames(cov_H1) <- parnames
  
  ### compute BF10 using Laplace approximation ###
  
  logml0 <- unname(.5 * log(2 * pi) +  .5 * log(var_H0) -
                     minl0(par = mode_H0, data = data, mu_beta = mu_beta,
                           sigma_beta = sigma_beta))
  logml1 <- unname(log(2 * pi) -
                     .5 * log(det_hessian_minl(par = mode_H1, data = data,
                                               sigma_beta = sigma_beta,
                                               sigma_psi = sigma_psi)) -
                     minl(par = mode_H1, data = data, mu_beta = mu_beta,
                          sigma_beta = sigma_beta, mu_psi = mu_psi,
                          sigma_psi = sigma_psi))
  logbf10 <- logml1 - logml0
  
  ### compute one-sided BFs ###
  
  method <- "log-is" # importance sampling with t-proposal based on
  # Laplace approximation to log transformed posterior
  
  r <- try({
    
    if (start_psi > 0) {
      start_xi_plus <- log(start_psi)
      start_xi_minus <- log(.1)
    } else if (start_psi < 0) {
      start_xi_plus <- log(.1)
      start_xi_minus <- log( - start_psi)
    } else {
      start_xi_plus <- 0
      start_xi_minus <- 0
    }
    
    # H+
    oplus <- stats::optim(par = c(start_beta, start_xi_plus), fn = minlplus,
                   gr = gradient_minlplus, data = data,
                   method = "BFGS", mu_beta = mu_beta,
                   sigma_beta = sigma_beta, mu_psi = mu_psi,
                   sigma_psi = sigma_psi)
    mode_Hplus <- oplus$par
    cov_Hplus <- inverse_hessian_minlplus(par = mode_Hplus,
                                          data = data,
                                          sigma_beta = sigma_beta,
                                          mu_psi = mu_psi,
                                          sigma_psi = sigma_psi)
    
    # proposal samples
    prop_samples_Hplus <- mvtnorm::rmvt(nsamples, delta = mode_Hplus,
                                        sigma = cov_Hplus, df = is_df)
    
    # importance weights
    weights_Hplus <- - apply_minlplus_cpp(prop_samples_Hplus, data$y1,
                                          data$y2, data$n1, data$n2,
                                          mu_beta, sigma_beta, mu_psi,
                                          sigma_psi) -
      mvtnorm::dmvt(prop_samples_Hplus, delta = mode_Hplus,
                    sigma = cov_Hplus, df = is_df, log = TRUE)
    logconst_Hplus <- stats::median(weights_Hplus)
    logmlplus <- log(mean(exp(weights_Hplus - logconst_Hplus))) +
      logconst_Hplus
    
    # compute BF+0
    logbfplus0 <- logmlplus - logml0
    
    # H-
    ominus <- stats::optim(par = c(start_beta, start_xi_minus), fn = minlminus,
                    gr = gradient_minlminus, data = data,
                    method = "BFGS", mu_beta = mu_beta,
                    sigma_beta = sigma_beta, mu_psi = mu_psi,
                    sigma_psi = sigma_psi)
    mode_Hminus <- ominus$par
    cov_Hminus <- inverse_hessian_minlminus(par = mode_Hminus,
                                            data = data,
                                            sigma_beta = sigma_beta,
                                            mu_psi = mu_psi,
                                            sigma_psi = sigma_psi)
    
    # proposal samples
    prop_samples_Hminus <- mvtnorm::rmvt(nsamples, delta = mode_Hminus,
                                         sigma = cov_Hminus, df = is_df)
    
    # importance weights
    weights_Hminus <- - apply_minlminus_cpp(prop_samples_Hminus,
                                            data$y1, data$y2,
                                            data$n1, data$n2, mu_beta,
                                            sigma_beta, mu_psi,
                                            sigma_psi) -
      mvtnorm::dmvt(prop_samples_Hminus, delta = mode_Hminus,
                    sigma = cov_Hminus, df = is_df, log = TRUE)
    logconst_Hminus <- stats::median(weights_Hminus)
    logmlminus <- log(mean(exp(weights_Hminus - logconst_Hminus))) +
      logconst_Hminus
    
    # compute BF-0
    logbfminus0 <- logmlminus - logml0
    
  }, silent = TRUE)
  
  if (inherits(r, "try-error") || any(is.na(c(logbfplus0, logbfminus0)))) {
    
    method <- "is-sn" # importance sampling to obtain unconstrained samples
    # then fit skew-normal and compute areas of interest
    
    mode_Hplus <- cov_Hplus <- NULL
    mode_Hminus <- cov_Hminus <- NULL
    
    # compute log of prior area larger than zero for psi
    logPriorPlus <- stats::pnorm(0, mean = mu_psi, sd = sigma_psi,
                          lower.tail = FALSE, log.p = TRUE)
    # compute log of prior area smaller than zero for psi
    logPriorMinus <- stats::pnorm(0, mean = mu_psi, sd = sigma_psi, log.p = TRUE)
    
    # obtain posterior samples via importance sampling
    
    # proposal samples
    prop_samples_H1 <- mvtnorm::rmvt(nsamples, delta = mode_H1,
                                     sigma = cov_H1, df = is_df)
    
    # importance weights
    weights_H1 <- - apply_minl_cpp(prop_samples_H1, data$y1, data$y2,
                                   data$n1, data$n2, mu_beta,
                                   sigma_beta, mu_psi, sigma_psi) -
      mvtnorm::dmvt(prop_samples_H1, delta = mode_H1, sigma = cov_H1,
                    df = is_df, log = TRUE)
    
    logconst_H1 <- stats::median(weights_H1)
    
    # normalized importance weights
    normalized_weights_H1 <- exp(weights_H1 - logconst_H1) /
      sum(exp(weights_H1 - logconst_H1))
    
    # resample according to normalize weights to get posterior samples
    index <- sample(1:nrow(prop_samples_H1), size = nsamples,
                    replace = TRUE, prob = normalized_weights_H1)
    post_samples_H1 <- prop_samples_H1[index,]
    colnames(post_samples_H1) <- parnames
    
    # fit skew-normal distribution and compute area smaller/larger zero
    psi <- post_samples_H1[,"psi"]
    fit <- sn::selm(formula = psi ~ 1, family = "SN",
                    data = data.frame(psi = psi))
    logPostMinus <- log(sn::psn(x = 0, xi = fit@param$dp[["xi"]],
                                omega = fit@param$dp[["omega"]],
                                alpha = fit@param$dp[["alpha"]]))
    logPostPlus <- log(1 - exp(logPostMinus))
    
    logbfplus0 <- logbf10 + logPostPlus - logPriorPlus
    logbfminus0 <- logbf10 + logPostMinus - logPriorMinus
    
  }
  
  out <- list(input = list(data = data,
                           prior_par = prior_par,
                           nsamples = nsamples,
                           is_df = is_df),
              method = method,
              bf = list(bf10 = exp(logbf10),
                        bfplus0 = exp(logbfplus0),
                        bfminus0 = exp(logbfminus0)))
  
  return(out)
  
  }
