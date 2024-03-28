#include "../inst/include/RIVSparse.h"

using namespace Rcpp;

// Base class
class BaseVCSC {
 public:
  virtual ~BaseVCSC() {}
};

// Templated derived class
template <typename T, typename U>
class DerivedVCSC : public BaseVCSC {
 public:
  IVSparse::VCSC<T, U> *mat;
};

// Your VCSC class
class VCSC {
 public:

  VCSC() {
    Rcout << "VCSC() called" << std::endl;
  }

  VCSC(const S4 &mat, std::string data_type, std::string index_type) {
    if (data_type == "int" && index_type == "int") {
      DerivedVCSC<int, int> *derived_vcsc_mat = new DerivedVCSC<int, int>();

      // initialize the VCSC matrix
      IntegerVector dims = mat.slot("Dim");
      IntegerVector i = mat.slot("i");
      IntegerVector p = mat.slot("p");
      IntegerVector x = mat.slot("x");

      int nrow = dims[0];
      int ncol = dims[1];
      int nnz = i.size();

      // get the data from the i, p, and x slots into a std::vector
      std::vector<int> i2(i.begin(), i.end());     // indices
      std::vector<int> p2(p.begin(), p.end());     // column pointers
      std::vector<int> x2(x.begin(), x.end());     // values

      // get pointers to the vectors
      int *i2_ptr = &i2[0];
      int *p2_ptr = &p2[0];
      int *x2_ptr = &x2[0];

      // initialize the VCSC matrix
      derived_vcsc_mat->mat = new IVSparse::VCSC<int, int>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);

      // Assign the derived_vcsc_mat to vcsc_mat
      vcsc_mat = derived_vcsc_mat;

      // print nnz using nonZeros() vcsc method
      int nnz_vcsc = derived_vcsc_mat->mat->nonZeros();
      Rcout << "nnz_vcsc: " << nnz_vcsc << std::endl;
    }
    // Add more conditions for other data types

  }

  VCSC(double x) {
    Rcout << "VCSC(double x) called" << std::endl;
  }

  ~VCSC() { delete vcsc_mat; }

  // Add methods here



 private:
  BaseVCSC *vcsc_mat;
};

// Expose the VCSC class
RCPP_MODULE(vcsc) {
  class_<VCSC>("VCSC")
    .constructor()
    .constructor<S4, std::string, std::string>()
    .constructor<double>()
    // Add methods here
    ;
}
