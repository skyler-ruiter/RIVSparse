#include "../inst/include/RIVSparse.h"

using namespace Rcpp;

// Namespace for package
namespace RIVSparse {

// VCSC class
template <typename T, typename U>
class VCSC {
  public:

  //* ---------------------Attributes--------------------- *//

  IVSparse::VCSC<T, U> *vcsc;
  bool columnMajor = true;

  //* ---------------------Constructors--------------------- *//
  //TODO: filename constructor and IVCSC constructor

  // default constructor
  VCSC() {}

  // constructor from a dg(C|R|T)Matrix (S4 object)
  VCSC(const S4 &mat) {
    IntegerVector i, p, Dim;
    NumericVector x;
    List Dimnames;

    // get sparse matrix type
    std::string mat_type = mat.attr("class");

    if (mat_type == "dgCMatrix") {
      i = mat.slot("i");
      p = mat.slot("p");
      x = mat.slot("x");
      Dim = mat.slot("Dim");
    } else if (mat_type == "dgRMatrix") {
      i = mat.slot("j");
      p = mat.slot("p");
      x = mat.slot("x");
      Dim = mat.slot("Dim");
      columnMajor = false;
    } else if (mat_type == "dgTMatrix") {
      i = mat.slot("i");
      p = mat.slot("j");
      x = mat.slot("x");
      Dim = mat.slot("Dim");
    } else {
      // print mat_type
      printf("Matrix type: %s\n", mat_type.c_str());

      stop("Invalid matrix type");
    }

    uint64_t nrow = Dim[0];
    uint64_t ncol = Dim[1];
    uint64_t nnz = x.size();

    // get pointers to the vectors
    int *i_ptr = &i[0];
    int *p_ptr = &p[0];
    double *x_ptr = &x[0];

    // make new poitners of type <T, U>
    U *i_ptr_new;
    U *p_ptr_new;
    T *x_ptr_new;

    // if type of <T, U> is the same as the original type, just use the original pointers
    if (std::is_same<T, double>::value && std::is_same<U, int>::value) {
      i_ptr_new = (U *) i_ptr;
      p_ptr_new = (U *) p_ptr;
      x_ptr_new = (T *) x_ptr;
    } else {
            // convert the pointers to the new type
      i_ptr_new = new U[nnz];
      p_ptr_new = new U[ncol + 1];
      x_ptr_new = new T[nnz];

      for (int k = 0; k < nnz; k++) {
        i_ptr_new[k] = i_ptr[k];
        x_ptr_new[k] = x_ptr[k];
      }

      for (int k = 0; k < ncol + 1; k++) {
        p_ptr_new[k] = p_ptr[k];
      }
    }

    if (mat_type == "dgTMatrix") {
      // construct a vector of tuples to store the triplet matrix
      std::vector<std::tuple<U, U, T>> entries;
      for (int k = 0; k < nnz; k++) {
        entries.push_back(std::make_tuple(i_ptr_new[k], p_ptr_new[k], x_ptr_new[k]));
      }

      // sort the entries by column
      std::sort(entries.begin(), entries.end(), [](const std::tuple<U, U, T> &a, const std::tuple<U, U, T> &b) {
        return std::get<1>(a) < std::get<1>(b);
      });

      vcsc = new IVSparse::VCSC<T, U>(entries, nrow, ncol, nnz);
      return;
    }

    vcsc = new IVSparse::VCSC<T, U>(x_ptr_new, i_ptr_new, p_ptr_new, nrow, ncol, nnz);
  }

  // destructor
  ~VCSC() { delete vcsc; }

  //* ---------------------Getters--------------------- *//
  //TODO: getters for specific VCSC data (vals, counts, indices, etc.)

  // coeff
  T coeff(const uint64_t i, const uint64_t j) { return vcsc->coeff(i, j); }

  // rows
  uint64_t rows() { return vcsc->rows(); }

  // cols
  uint64_t cols() { return vcsc->cols(); }

  // inner dimension
  uint64_t innerDim() { return vcsc->innerSize(); }

  // outer dimension
  uint64_t outerDim() { return vcsc->outerSize(); }

  // nnz
  uint64_t nnz() { return vcsc->nonZeros(); }

  // byte size
  uint64_t byteSize() { return vcsc->byteSize(); }

  // get the column major attribute
  bool getColumnMajor() { return vcsc->isColumnMajor(); }

  //* ---------------------Converters--------------------- *//
  //TODO: converters to dg(C|R|T)Matrix and IVCSC


  //* ---------------------Calculations--------------------- *//
  //TODO: col/rowSum, max/min, trace, sum, norm, etc.


  //* ---------------------Utility--------------------- *//
  //TODO: file I/O and cleaner print method

  // R print method
  void print() {
    // if more than 100 cols or rows print summary and first 10 rows and cols
    if (vcsc->rows() > 100 || vcsc->cols() > 100) {
      Rcout << "VCSC matrix with " << vcsc->rows() << " rows and " << vcsc->cols() << " columns" << std::endl;
      Rcout << "First 10 rows and columns:" << std::endl;
      for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
          Rcout << vcsc->coeff(i, j) << " ";
        }
        Rcout << std::endl;
      }
    } else {
      Rcout << "VCSC matrix with " << vcsc->rows() << " rows and " << vcsc->cols() << " columns" << std::endl;
      for (int i = 0; i < vcsc->rows(); i++) {
        for (int j = 0; j < vcsc->cols(); j++) {
          Rcout << vcsc->coeff(i, j) << " ";
        }
        Rcout << std::endl;
      }
    }
  }

  //* ---------------------Matrix Manipulation--------------------- *//
  //TODO: transpose (and in place transpose)

  // append a matrix to the right of the current matrix
  void append(SEXP vcsc_mat) {
    // get the vcsc object
    XPtr<IVSparse::VCSC<T, U>> vcsc_xptr(vcsc_mat);
    IVSparse::VCSC<T, U> *vcsc_obj = vcsc_xptr;

    // append the vcsc object
    vcsc->append(*vcsc_obj);
  }

  // get a slice of the matrix columns/rows
  SEXP slice(uint64_t start, uint64_t end) {
    // create a new XPtr object
    XPtr<IVSparse::VCSC<T, U>> slice_vcsc_ptr(new IVSparse::VCSC<T, U>(vcsc->slice(start, end)), true);
    return slice_vcsc_ptr;
  }

  //* --------------------Matrix Operations (BLAS)-------------------- *//
  //TODO: return scalar, SpMM, and SpMV

  // scale the matrix by a scalar in place
  void scale(const T scalar) { vcsc->operator*=(scalar); }

  //* ---------------------Operator Overloads?--------------------- *//
  //TODO: equality/inequality, assignment, etc.
  
  //* ---------------------R6 Methods--------------------- *//

  // make method to return underlying vcsc object as a smart pointer
  SEXP get_vcsc_xptr() {
    XPtr<IVSparse::VCSC<T, U>> vcsc_xptr(vcsc, true);
    return vcsc_xptr;
  }

  void set_vcsc_ptr(SEXP vcsc_obj) {
    XPtr<IVSparse::VCSC<T, U>> vcsc_xptr(vcsc_obj);
    vcsc = new IVSparse::VCSC<T, U>(*vcsc_xptr);
  }

  //* ---------------------Test Methods--------------------- *//
  //TODO: clean up test methods

  //! --Successs-- !// taking in a vcsc object as argument!
  // void test2(Environment &vcsc) {
  //   Rprintf("Printing from test2\n");
  //   Function print = vcsc["print"];
  //   print();
  // }

  // test method that returns an R6 VCSC object
  //! ends up as another shallow copy of the original object
  // Environment test_return(Environment &vcsc_env) {
  //   // create a new environment
  //   Environment package_env("package:RIVSparse");
  //   Environment class_env = package_env["VCSC"];
  //   Function new_vcsc = class_env["new"];

  //   // create a new vcsc object providing the vcsc object as an argument
  //   Environment new_vcsc_env = new_vcsc(vcsc_env);
  //   return new_vcsc_env;
  // }

};
}

// Define the types for each data and index type
//TODO: add the rest of the data types
typedef RIVSparse::VCSC<int, int> VCSC_INT_INT;               // int    int
typedef RIVSparse::VCSC<int, uint64_t> VCSC_INT_UINT64;       // int    uint64_t
typedef RIVSparse::VCSC<double, int> VCSC_DOUBLE_INT;         // double int
typedef RIVSparse::VCSC<double, uint64_t> VCSC_DOUBLE_UINT64; // double uint64_t

//* int int *//
RCPP_MODULE(vcsc_int_int) {
  class_<VCSC_INT_INT>("VCSC_INT_INT")
  // constructors
  .constructor()
  .constructor<S4>()
  // getters
  .method("coeff", &VCSC_INT_INT::coeff)
  // converters
  // calculations
  // utility
  .method("print", &VCSC_INT_INT::print)
  // matrix manipulation
  .method("slice", &VCSC_INT_INT::slice)
  .method("append", &VCSC_INT_INT::append)
  // matrix operations
  .method("scale", &VCSC_INT_INT::scale)
  // operator overloads
  // R6 methods
  .method("get_vcsc_xptr", &VCSC_INT_INT::get_vcsc_xptr)
  .method("set_vcsc_ptr", &VCSC_INT_INT::set_vcsc_ptr)
  ;
}

//* int uint64_t *//
RCPP_MODULE(vcsc_int_uint64) {
  class_<VCSC_INT_UINT64>("VCSC_INT_UINT64")
  // constructors
  .constructor()
  .constructor<S4>()
  // getters
  .method("coeff", &VCSC_INT_UINT64::coeff)
  // converters
  // calculations
  // utility
  .method("print", &VCSC_INT_UINT64::print)
  // matrix manipulation
  .method("append", &VCSC_INT_UINT64::append)
  .method("slice", &VCSC_INT_UINT64::slice)
  // matrix operations
  .method("scale", &VCSC_INT_UINT64::scale)
  // operator overloads
  // R6 methods
  .method("get_vcsc_xptr", &VCSC_INT_UINT64::get_vcsc_xptr)
  .method("set_vcsc_ptr", &VCSC_INT_UINT64::set_vcsc_ptr)
  ;
}

//* double int *//
RCPP_MODULE(vcsc_double_int) {
  class_<VCSC_DOUBLE_INT>("VCSC_DOUBLE_INT")
  // constructors
  .constructor()
  .constructor<S4>()
  // getters
  .method("coeff", &VCSC_DOUBLE_INT::coeff)
  // converters
  // calculations
  // utility
  .method("print", &VCSC_DOUBLE_INT::print)
  // matrix manipulation
  .method("append", &VCSC_DOUBLE_INT::append)
  .method("slice", &VCSC_DOUBLE_INT::slice)
  // matrix operations
  .method("scale", &VCSC_DOUBLE_INT::scale)
  // operator overloads
  // R6 methods
  .method("get_vcsc_xptr", &VCSC_DOUBLE_INT::get_vcsc_xptr)
  .method("set_vcsc_ptr", &VCSC_DOUBLE_INT::set_vcsc_ptr)
  ;
}

//* double uint64_t *//
RCPP_MODULE(vcsc_double_uint64) {
  class_<VCSC_DOUBLE_UINT64>("VCSC_DOUBLE_UINT64")
  // constructors
  .constructor()
  .constructor<S4>()
  // getters
  .method("coeff", &VCSC_DOUBLE_UINT64::coeff)
  // converters
  // calculations
  // utility
  .method("print", &VCSC_DOUBLE_UINT64::print)
  // matrix manipulation
  .method("append", &VCSC_DOUBLE_UINT64::append)
  .method("slice", &VCSC_DOUBLE_UINT64::slice)
  // matrix operations
  .method("scale", &VCSC_DOUBLE_UINT64::scale)
  // operator overloads
  // R6 methods
  .method("get_vcsc_xptr", &VCSC_DOUBLE_UINT64::get_vcsc_xptr)
  .method("set_vcsc_ptr", &VCSC_DOUBLE_UINT64::set_vcsc_ptr)
  ;
}