#include <Rcpp.h>

#include <iostream>

// enable openmp
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//' Compute ACF using OpenMP
//'
//' ACF is computed in 'cov'-type.
//' @param x The input data.
//' @param tau_max The max correlation lag.
// [[Rcpp::export]]
NumericVector fastacf(NumericVector x, unsigned int tau_max) {
  unsigned int n = x.length();
  // compute periodic mean
  double mu = 0;
  for (unsigned int i=0; i < n; ++i) {
    mu += x[i];
  }
  mu /= n;
  // shift coordinates using periodic boundary corrections
  for (unsigned int i=0; i < n; ++i) {
    x[i] = x[i] - mu;
  }
  // compute ACF
  unsigned int tau;
  unsigned int i;
  std::vector<double> acf(tau_max, 0);
  #pragma omp parallel for default(none)\
    private(tau,i)\
    firstprivate(tau_max,n)\
    shared(acf,x)\
    schedule(static,16)
  for (unsigned int tau=0; tau < tau_max; ++tau) {
    for (unsigned int i=0; i < n-tau; ++i) {
      acf[tau] += x[i+tau] * x[i];
    }
    acf[tau] /= (n-tau);
  }
  return NumericVector(acf.begin(), acf.end());
}
