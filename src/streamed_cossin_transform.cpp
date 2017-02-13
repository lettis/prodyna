#include <Rcpp.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iostream>

using namespace Rcpp;

//' computes a cos/sin transform of angles on the fly
//'
//' @param fname_in Input filename with angles.
//' @param fname_out Output filename with cos/sin transformed values.
//' @param is_deg Angles are expressed in degrees. default: true.
// [[Rcpp::export]]
void
streamed_cossin_transform(std::string fname_in
                        , std::string fname_out
                        , bool is_deg = true) {
  std::ifstream ifs(fname_in);
  std::ofstream ofs(fname_out);
  while (ifs.good() && ofs.good()) {
    std::string buf;
    std::getline(ifs, buf);
    if ( ! buf.empty()) {
      std::stringstream ss(buf);
      while (ss.good()) {
        float value;
        ss >> value;
        if (! ss.fail()) {
          if (is_deg) {
            value *= (M_PI / 180.0);
          }
          ofs << " " << std::cos(value) << " " << std::sin(value);
        }
      }
      ofs << "\n";
    }
  }
}
