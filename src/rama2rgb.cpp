#include <Rcpp.h>

#include <iostream>

// enable openmp
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

float
pow2(float x) {
  return x*x;
}

//' Compute RGB color code from sampling of Ramachandran plot
//'
//' phi and psi values are defined on range [-180, 180].
//' @param phis Sampling of phi values.
//' @param psis Sampling of psi values.
// [[Rcpp::export]]
NumericVector rama2rgb(NumericVector phis
                     , NumericVector psis) {
  unsigned int n = phis.length();
  std::vector<float> v_phis(phis.begin()
                          , phis.end());
  std::vector<float> v_psis(psis.begin()
                          , psis.end());

  // color weight for given refernce value
  auto weight = [&](unsigned int i
                  , float ref_phi
                  , float ref_psi) -> float {
    constexpr float scale = std::sqrt(2) * 180.0;
    float value = std::min(pow2(v_phis[i]-ref_phi)
                         , std::min(pow2(v_phis[i]-ref_phi+360)
                                  , pow2(v_phis[i]-ref_phi-360)))
                + std::min(pow2(v_psis[i]-ref_psi)
                         , std::min(pow2(v_psis[i]-ref_psi+360)
                                  , pow2(v_psis[i]-ref_psi-360)));
    value = 1 - (std::sqrt(value) / scale);
    return value;
  };


  unsigned int i;
  float r=0;
  float g=0;
  float b=0;
  #pragma omp parallel for\
    default(none)\
    private(i)\
    firstprivate(n,weight)\
    reduction(+:r,g,b)
  for (i=0; i < n; ++i) {
    r += weight(i, -120, 120);
    g += weight(i, -60, -60);
    b += weight(i, 60, 60);
  }
  std::vector<float> rgb{r,g,b};
  for (float& clr: rgb) {
    clr /= n;
  }
  return NumericVector(rgb.begin(), rgb.end());
}
