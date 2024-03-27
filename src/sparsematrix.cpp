#include "../inst/include/RIVSparse.h"

using namespace Rcpp;

class Uniform {
 public:
  Uniform(double min_, double max_) : min(min_), max(max_) {}
  NumericVector draw(int n) const {
    RNGScope scope;
    return runif(n, min, max);
  }
  double min, max;
};

double uniformRange(Uniform* w) { return w->max - w->min; }

RCPP_MODULE(unif_module) {
  class_<Uniform>("Uniform")

      .constructor<double, double>()

      .field("min", &Uniform::min)
      .field("max", &Uniform::max)

      .method("draw", &Uniform::draw)
      .method("range", &uniformRange);
}





// [[Rcpp::export]]
void convertSparse(S4 mat) {

  IntegerVector dims = mat.slot("Dim");
  IntegerVector i = mat.slot("i");
  IntegerVector p = mat.slot("p");
  NumericVector x = mat.slot("x");

  int nrow = dims[0];
  int ncol = dims[1];

  // get the data from the i, p, and x slots into a std::vector
  std::vector<int> i2(i.begin(), i.end());
  std::vector<int> p2(p.begin(), p.end());
  std::vector<double> x2(x.begin(), x.end());

  // print the data
  Rcout << "nrow: " << nrow << std::endl;
  Rcout << "ncol: " << ncol << std::endl;

  Rcout << "i: ";
  for (int j = 0; j < i2.size(); j++) {
    Rcout << i2[j] << " ";
  }
}
