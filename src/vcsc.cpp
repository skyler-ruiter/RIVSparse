#include "../inst/include/RIVSparse.h"

using namespace Rcpp;
namespace RIVSparse {

template <typename T, typename U>
class VCSC {
  public:

  // Attributes
  IVSparse::VCSC<T, U> *vcsc;
  bool columnMajor = true;

  // Methods and Constructors
  VCSC() {}

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

  ~VCSC() { delete vcsc; }

  // coeff
  T coeff(const uint64_t i, const uint64_t j) {
    return vcsc->coeff(i, j);
  }

  // scale the matrix by a scalar returning a new matrix
  void scale(const T scalar) {
    vcsc->operator*=(scalar);
  }

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


  // slice method returning a new VCSC object
  VCSC<T, U> slice(const uint64_t start, const uint64_t end) {
    return VCSC<T, U>(vcsc->slice(start, end));
  }

  //* ---------------------Test Methods--------------------- *//

  // test method for template specializations

  // test method for returning a vcsc object using an environment object
  // Environment test(double scalar) {
  //   // create a new environment
  //   Environment package_env("package:RIVSparse");
  //   Environment class_env = package_env["VCSC"];
  //   Function new_vcsc = class_env["new"];

  //   // create a new vcsc object
  //   Environment vcsc_env = new_vcsc();
  // }

  //! --Successs-- !// taking in a vcsc object as argument!
  void test2(Environment &vcsc) {
    Rprintf("Printing from test2\n");
    Function print = vcsc["print"];
    print();
  }


};
}

typedef RIVSparse::VCSC<int, int> VCSC_INT_INT;
typedef RIVSparse::VCSC<int, uint64_t> VCSC_INT_UINT64;
typedef RIVSparse::VCSC<double, int> VCSC_DOUBLE_INT;
typedef RIVSparse::VCSC<double, uint64_t> VCSC_DOUBLE_UINT64;


// *template specializations for Rcpp

// template <>
// VCSC_INT_INT VCSC_INT_INT::test() {
//   return VCSC_INT_INT();
// }

// template <>
// VCSC_INT_UINT64 VCSC_INT_UINT64::test() {
//   return VCSC_INT_UINT64();
// }

// template <>
// VCSC_DOUBLE_INT VCSC_DOUBLE_INT::test() {
//   return VCSC_DOUBLE_INT();
// }

// template <>
// VCSC_DOUBLE_UINT64 VCSC_DOUBLE_UINT64::test() {
//   return VCSC_DOUBLE_UINT64();
// }





RCPP_MODULE(vcsc_int_int) {
  class_<VCSC_INT_INT>("VCSC_INT_INT")
  .constructor()
  .constructor<S4>()
  .method("coeff", &VCSC_INT_INT::coeff)
  .method("scale", &VCSC_INT_INT::scale)
  .method("print", &VCSC_INT_INT::print)
  // .method("test", &VCSC_INT_INT::test)
  .method("test2", &VCSC_INT_INT::test2)
  ;
}

RCPP_MODULE(vcsc_int_uint64) {
  class_<VCSC_INT_UINT64>("VCSC_INT_UINT64")
  .constructor()
  .constructor<S4>()
  .method("coeff", &VCSC_INT_UINT64::coeff)
  .method("scale", &VCSC_INT_UINT64::scale)
  .method("print", &VCSC_INT_UINT64::print)
  // .method("test", &VCSC_INT_UINT64::test)
  .method("test2", &VCSC_INT_UINT64::test2)
  ;
}

RCPP_MODULE(vcsc_double_int) {
  class_<VCSC_DOUBLE_INT>("VCSC_DOUBLE_INT")
  .constructor()
  .constructor<S4>()
  .method("coeff", &VCSC_DOUBLE_INT::coeff)
  .method("scale", &VCSC_DOUBLE_INT::scale)
  .method("print", &VCSC_DOUBLE_INT::print)
  // .method("test", &VCSC_DOUBLE_INT::test)
  .method("test2", &VCSC_DOUBLE_INT::test2)
  ;
}

RCPP_MODULE(vcsc_double_uint64) {
  class_<VCSC_DOUBLE_UINT64>("VCSC_DOUBLE_UINT64")
  .constructor()
  .constructor<S4>()
  .method("coeff", &VCSC_DOUBLE_UINT64::coeff)
  .method("scale", &VCSC_DOUBLE_UINT64::scale)
  .method("print", &VCSC_DOUBLE_UINT64::print)
  // .method("test", &VCSC_DOUBLE_UINT64::test)
  .method("test2", &VCSC_DOUBLE_UINT64::test2)
  ;
}