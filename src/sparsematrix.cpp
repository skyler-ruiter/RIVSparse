#include "../inst/include/RIVSparse.h"

using namespace Rcpp;


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

// make a function that takes in data from the R6 object wrapper and returns an XPtr to the C++ object
// [[Rcpp::export]]
Rcpp::XPtr<IVSparse::VCSC<double, int>> create_vcsc(S4 mat) {
  IntegerVector dims = mat.slot("Dim");
  IntegerVector i = mat.slot("i");
  IntegerVector p = mat.slot("p");
  NumericVector x = mat.slot("x");

  int nrow = dims[0];
  int ncol = dims[1];

  // get pointers to the data
  double* data_ptr = x.begin();
  int* ind_ptr = i.begin();
  int* ptr_ptr = p.begin();

  size_t nnz = x.size();

  IVSparse::VCSC<double, int>* ptr = new IVSparse::VCSC<double, int>(data_ptr, ind_ptr, ptr_ptr, nrow, ncol, nnz);

  return Rcpp::XPtr<IVSparse::VCSC<double, int>>(ptr);
}

RcppExport SEXP RcppExports(Rcpp::Environment env) {
  env.function("create_vcsc", &create_vcsc);
  return env;
}