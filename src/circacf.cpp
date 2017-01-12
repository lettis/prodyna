#include <Rcpp.h>
using namespace Rcpp;

//' Compute ACF of circular data
//'
//' Computes the autocorrelation of circular data, defined on
//' the range [0, 2pi], expressed in radians.
//' ACF is computed in 'cov'-type.
//' @param x The input data.
//' @param tau_max The max correlation lag.
// [[Rcpp::export]]
NumericVector circacf(NumericVector x, unsigned int tau_max) {
  NumericVector rho(tau_max);

  //TODO implement for radians and [0,2pi] range

  return rho;
}
