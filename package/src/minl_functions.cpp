// load Rcpp
#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// [[Rcpp::export]]
double minl_cpp(double beta, double psi, int y1, int y2, int n1,
                int n2, double mu_beta, double sigma_beta,
                double mu_psi, double sigma_psi) {
  
  
  // log probabilities
  double logp1 = beta - 0.5 * psi - R::log1pexp(beta - 0.5 * psi);
  double log1minp1 = - log1pexp(beta - 0.5 * psi);
  double logp2 = beta + 0.5 * psi - log1pexp(beta + 0.5 * psi);
  double log1minp2 = - log1pexp(beta + 0.5 * psi);
  
  double out = y1 * logp1 + (n1 - y1) * log1minp1 +
    y2 * logp2 + (n2 - y2) * log1minp2 +
    R::dnorm(beta, mu_beta, sigma_beta, 1) +
    R::dnorm(psi, mu_psi, sigma_psi, 1);
  
  return(- out);
  
}

// [[Rcpp::export]]
NumericVector apply_minl_cpp(NumericMatrix s, int y1, int y2,
                            int n1, int n2, double mu_beta,
                            double sigma_beta, double mu_psi,
                            double sigma_psi) {
  
  int nsamples = s.nrow();
  NumericVector out(nsamples);
  
  for (int i = 0; i < nsamples; i++) {
    
    out(i) = minl_cpp(s(i,0), s(i,1), y1, y2, n1, n2, mu_beta,
        sigma_beta, mu_psi, sigma_psi);
    
  }
  
  return out;
  
}

// [[Rcpp::export]]
double minlminus_cpp(double beta, double xi, int y1, int y2, int n1,
                     int n2, double mu_beta, double sigma_beta,
                     double mu_psi, double sigma_psi) {
  
  double psi = - exp(xi);
  
  // log probabilities
  double logp1 = beta - 0.5 * psi - R::log1pexp(beta - 0.5 * psi);
  double log1minp1 = - log1pexp(beta - 0.5 * psi);
  double logp2 = beta + 0.5 * psi - log1pexp(beta + 0.5 * psi);
  double log1minp2 = - log1pexp(beta + 0.5 * psi);
  
  double out = y1 * logp1 + (n1 - y1) * log1minp1 +
    y2 * logp2 + (n2 - y2) * log1minp2 +
    R::dnorm(beta, mu_beta, sigma_beta, 1) +
    R::dnorm(psi, mu_psi, sigma_psi, 1) -
    R::pnorm(0, mu_psi, sigma_psi, 1, 1) + xi;
  
  return(- out);
  
}

// [[Rcpp::export]]
NumericVector apply_minlminus_cpp(NumericMatrix s, int y1, int y2,
                                  int n1, int n2, double mu_beta,
                                  double sigma_beta, double mu_psi,
                                  double sigma_psi) {
  
  int nsamples = s.nrow();
  NumericVector out(nsamples);
  
  for (int i = 0; i < nsamples; i++) {
   
   out(i) = minlminus_cpp(s(i,0), s(i,1), y1, y2, n1, n2, mu_beta,
       sigma_beta, mu_psi, sigma_psi);
    
  }
  
  return out;
  
}

// [[Rcpp::export]]
double minlplus_cpp(double beta, double xi, int y1, int y2, int n1,
                    int n2, double mu_beta, double sigma_beta,
                    double mu_psi, double sigma_psi) {
  
  double psi = exp(xi);
  
  // log probabilities
  double logp1 = beta - 0.5 * psi - R::log1pexp(beta - 0.5 * psi);
  double log1minp1 = - log1pexp(beta - 0.5 * psi);
  double logp2 = beta + 0.5 * psi - log1pexp(beta + 0.5 * psi);
  double log1minp2 = - log1pexp(beta + 0.5 * psi);
  
  double out = y1 * logp1 + (n1 - y1) * log1minp1 +
    y2 * logp2 + (n2 - y2) * log1minp2 +
    R::dnorm(beta, mu_beta, sigma_beta, 1) +
    R::dnorm(psi, mu_psi, sigma_psi, 1) -
    R::pnorm(0, mu_psi, sigma_psi, 0, 1) + xi;
  
  return(- out);
  
}

// [[Rcpp::export]]
NumericVector apply_minlplus_cpp(NumericMatrix s, int y1, int y2,
                                 int n1, int n2, double mu_beta,
                                 double sigma_beta, double mu_psi,
                                 double sigma_psi) {
  
  int nsamples = s.nrow();
  NumericVector out(nsamples);
  
  for (int i = 0; i < nsamples; i++) {
    
    out(i) = minlplus_cpp(s(i,0), s(i,1), y1, y2, n1, n2, mu_beta,
        sigma_beta, mu_psi, sigma_psi);
    
  }
  
  return out;
  
}
