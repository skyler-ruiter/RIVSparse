#include "../inst/include/RIVSparse.h"

using namespace Rcpp;

std::string hello() {
  return "hello";
}

void bla() {
  Rprintf("hello\\n");
}

void bla2( int x, double y) {
  Rprintf("hello (x = %d, y = %5.2f)\\n", x, y);
}

class World {
  public:
    World() : msg("hello") {}
    void set(std::string msg) { this->msg = msg; }
    std::string greet() { return msg; }

  private:
    std::string msg;
};

RCPP_MODULE(yada) {
  using namespace Rcpp;
  function("hello" , &hello);
  function("bla" , &bla);
  function("bla2" , &bla2);
  class_<World>("World")
    .constructor()
    .method("greet", &World::greet)
    .method("set", &World::set)
  ;
}




// [[Rcpp::export]]
void convertSparse(S4 &mat) {

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
