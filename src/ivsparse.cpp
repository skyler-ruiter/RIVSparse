#include "../inst/include/RIVSparse.h"

using namespace Rcpp;

// Base class
class BaseVCSC {
 public:
  virtual ~BaseVCSC() {}
  virtual double coeff(int i, int j) = 0;
  virtual NumericMatrix mult(NumericMatrix dense_mat) = 0;
};

// Templated derived class
template <typename T, typename U>
class DerivedVCSC : public BaseVCSC {
 public:
  IVSparse::VCSC<T, U> *mat;

  double coeff(int i, int j) {
    return mat->coeff(i, j);
  }

  NumericMatrix mult(NumericMatrix dense_mat) {
    // make a new dense matrix to store the result
    NumericMatrix result(dense_mat.ncol(), dense_mat.nrow());

    // transpose mat
    NumericMatrix mat_t = transpose(dense_mat);

    int outerDim = mat->cols();

    for (uint32_t i = 0; i < outerDim; ++i) {
      for (typename IVSparse::VCSC<T, U>::InnerIterator it(*mat, i); it; ++it) {
        result.column(it.getIndex()) = result.column(it.getIndex()) + (mat_t.column(i) * it.value());
      }
    }

    // return result transposed
    return transpose(result);
  }
};

// Your VCSC class
class VCSC {
 public:

  VCSC() {
    Rcout << "VCSC() called" << std::endl;
  }

  VCSC(const S4 &mat) {
    // create the vectors
    IntegerVector i, p, Dim;
    NumericVector x;
    List Dimnames;

    // check if the matrix is a dgCMatrix
    if (mat.inherits("dgCMatrix")) {
      i = mat.slot("i");
      p = mat.slot("p");
      x = mat.slot("x");
      Dim = mat.slot("Dim");
    } else {
      stop("Invalid matrix type");
    }

    uint64_t nrow = Dim[0];
    uint64_t ncol = Dim[1];
    uint64_t nnz = i.size();

    // get pointers to the vectors
    int *i_ptr = &i[0];
    int *p_ptr = &p[0];
    double *x_ptr = &x[0];

    // initialize the VCSC matrix
    DerivedVCSC<double, int> *derived_vcsc_mat = new DerivedVCSC<double, int>();

    derived_vcsc_mat->mat = new IVSparse::VCSC<double, int>(x_ptr, i_ptr, p_ptr, nrow, ncol, nnz);

    // init vcsc_mat
    vcsc_mat = derived_vcsc_mat;
    
    // print the matrix
    for (int i = 0; i < ncol; i++) {
      for (int j = 0; j < nrow; j++) {
        Rcout << derived_vcsc_mat->mat->coeff(i, j) << " ";
      }
      Rcout << std::endl;
    }
  }

  VCSC(const S4 &mat, std::string data_type, std::string index_type) {

    // create the vectors
    IntegerVector i, p, Dim;
    NumericVector x;
    List Dimnames;

    // check if the matrix is a dgCMatrix
    if (mat.inherits("dgCMatrix")) {
      i = mat.slot("i");
      p = mat.slot("p");
      x = mat.slot("x");
      Dim = mat.slot("Dim");
    } else {
      stop("Invalid matrix type");
    }

    uint64_t nrow = Dim[0];
    uint64_t ncol = Dim[1];
    uint64_t nnz = i.size();

    // get pointers to the vectors
    int *i_ptr = &i[0];
    int *p_ptr = &p[0];
    double *x_ptr = &x[0];

    if (data_type == "int" && index_type == "int") {
      DerivedVCSC<int, int> *derived_vcsc_mat = new DerivedVCSC<int, int>();

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

      // init vcsc_mat
      vcsc_mat = derived_vcsc_mat;

      // print nnz using nonZeros() vcsc method
      int nnz_vcsc = derived_vcsc_mat->mat->nonZeros();
      Rcout << "nnz_vcsc: " << nnz_vcsc << std::endl;

      // iterate through the vcsc matrix and print each vlaue
      for (int i = 0; i < ncol; i++) {
        for (int j = 0; j < nrow; j++) {
          Rcout << derived_vcsc_mat->mat->coeff(i, j) << " ";
        }
        Rcout << std::endl;
      }
    } 
    // else if (data_type == "" && index_type == "") {
    //   /* code */
    // } else if (data_type == "" && index_type == "") {
    //   /* code */
    // } else if (data_type == "" && index_type == "") {
    //   /* code */
    // } else {
    //   stop("Invalid data_type or index_type");
    // }
  }

  ~VCSC() { delete vcsc_mat; }

  // Add methods here
  double coeff(int i, int j) {
    return vcsc_mat->coeff(i, j);
  }

  NumericMatrix mult(NumericMatrix dense_mat) {
    return vcsc_mat->mult(dense_mat);
  }


 private:
  BaseVCSC *vcsc_mat;
};

// Expose the VCSC class
RCPP_MODULE(vcsc) {
  class_<VCSC>("VCSC")
    .constructor()
    .constructor<S4>()
    .constructor<S4, std::string, std::string>()
    // Add methods here
    .method("coeff", &VCSC::coeff)
    .method("mult", &VCSC::mult)
    ;
}
