#include "../inst/include/RIVSparse.h"

using namespace Rcpp;

class VCSC;

//! Base class
class BaseVCSC {
 public:

  // Constructor
  virtual ~BaseVCSC() {}

  //* Virtual methods
  /****************************************************
   *                                                   *
   *                                                   *
   *                   Converters                      *
   *                                                   *
   *                                                   *
   *****************************************************/
  // virtual S4 to_dgCMatrix() = 0;
  // virtual S4 to_dgRMatrix() = 0;
  // virtual S4 to_dgTMatrix() = 0;
  // virtual IVCSC to_ivcsc() = 0;

  /****************************************************
   *                                                   *
   *                                                   *
   *                      Getters                      *
   *                                                   *
   *                                                   *
   *****************************************************/
  virtual double coeff(int i, int j) = 0;
  virtual int rows() = 0;
  virtual int cols() = 0;
  virtual int nnz() = 0;
  virtual int innerdim() = 0;
  virtual int outerdim() = 0;
  virtual int bytesize() = 0;


  /****************************************************
   *                                                   *
   *                                                   *
   *               Matrix Manipulation                 *
   *                                                   *
   *                                                   *
   *****************************************************/
  virtual void transpose() = 0;
  // virtual VCSC slice(int start, int end) = 0;
  virtual void append(VCSC& other) = 0;
  virtual NumericMatrix mult(NumericMatrix dense_mat) = 0;

  /*****************************************************
   *                                                   *
   *                                                   *
   *                     Operators                     *
   *                                                   *
   *                                                   *
   *****************************************************/
  // virtual VCSC operator*(double scalar) = 0;
  // virtual VCSC operator*(NumericMatrix dense_mat) = 0;
  // virtual bool operator==(VCSC vcsc) = 0;
  // virtual bool operator!=(VCSC vcsc) = 0;
  // virtual VCSC operator=(VCSC vcsc) = 0;

  /*****************************************************
   *                                                   *
   *                                                   *
   *                       Misc                        *
   *                                                   *
   *                                                   *
   *****************************************************/
  // virtual void print() = 0;
  // virtual void save(std::string filename) = 0;
  // virtual VCSC load(std::string filename) = 0;
};






//! Templated derived class
template <typename T, typename U, bool columnMajor = 1>
class DerivedVCSC : public BaseVCSC {
 
 public:
  // Underlying VCSC matrix
  IVSparse::VCSC<T, U, columnMajor> *mat; // pointer to the VCSC matrix

  /****************************************************
   *                                                   *
   *                                                   *
   *                   Converters                      *
   *                                                   *
   *                                                   *
   *****************************************************/
  // S4 to_dgCMatrix() {}
  // S4 to_dgRMatrix() {}
  // S4 to_dgTMatrix() {}
  // IVCSC to_ivcsc() {}

  /******************************************************
   *                                                    *
   *                                                    *
   *                      Getters                       *
   *                                                    *
   *                                                    *
   *****************************************************/
  double coeff(int i, int j) { return mat->coeff(i, j); }
  int rows() { return mat->rows(); }
  int cols() { return mat->cols(); }
  int nnz() { return mat->nonZeros(); }
  int innerdim() { return mat->innerSize(); }
  int outerdim() { return mat->outerSize(); }
  int bytesize() { return mat->byteSize(); }

  /****************************************************
   *                                                   *
   *                                                   *
   *               Matrix Manipulation                 *
   *                                                   *
   *                                                   *
   *****************************************************/
  void transpose() {
    // transpose the matrix
    mat->inPlaceTranspose();
  }
  
  // VCSC slice(int start, int end) {}

  void append(VCSC& other) override;

  NumericMatrix mult(NumericMatrix dense_mat) {
    // make a new dense matrix to store the result
    NumericMatrix result(dense_mat.ncol(), dense_mat.nrow());

    // transpose mat
    NumericMatrix mat_t = Rcpp::transpose(dense_mat);

    std::vector<std::mutex> mutexList(mat->rows());

    int outerDim = mat->cols();

    #pragma omp parallel for
    for (uint32_t i = 0; i < outerDim; ++i) {
      for (typename IVSparse::VCSC<T, U, columnMajor>::InnerIterator it(*mat, i); it; ++it) {
        std::lock_guard<std::mutex> lock(mutexList[it.getIndex()]);
        result.column(it.getIndex()) = result.column(it.getIndex()) + (mat_t.column(i) * it.value());
      }
    }
    // return result transposed
    return Rcpp::transpose(result);
  }

  /*****************************************************
   *                                                   *
   *                                                   *
   *                     Operators                     *
   *                                                   *
   *                                                   *
   *****************************************************/
  // VCSC operator*(double scalar) {}
  // VCSC operator*(NumericMatrix dense_mat) {}
  // bool operator==(VCSC vcsc) {}
  // bool operator!=(VCSC vcsc) {}
  // VCSC operator=(VCSC vcsc) {}

  /*****************************************************
   *                                                   *
   *                                                   *
   *                       Misc                        *
   *                                                   *
   *                                                   *
   *****************************************************/
  // void print() {}
  // void save(std::string filename) {}
  // VCSC load(std::string filename) {}

  // IVSparse::VCSC<T, U, columnMajor> *getMat() override {
  //   return mat;
  // }
};




//! Your VCSC class
class VCSC {
 public:

 /****************************************************
 *                                                   *
 *                                                   *
 *                   Constructors                    *
 *                                                   *
 *                                                   *
 *****************************************************/

  // Default constructor
  VCSC() {}

  //! dgTMatrix->VCSC conversion is not great, needs to copy because not sorted COO
  // Constructor to convert a dg(C/R/T)Matrix to VCSC
  VCSC(const S4 &mat) {  // take a dg(C/R/T)Matrix and convert it to VCSC<int, double>
    // create the vectors
    IntegerVector i, p, Dim;
    NumericVector x;
    List Dimnames;
    bool columnMajor = true;

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

    if (mat_type == "dgTMatrix") {
      DerivedVCSC<double, int> *derived_vcsc_mat = new DerivedVCSC<double, int>();
      // construct a vector of tuples to store the triplet matrix
      std::vector<std::tuple<int, int, double>> entries;
      for (int k = 0; k < nnz; k++) {
        entries.push_back(std::make_tuple(i_ptr[k], p_ptr[k], x_ptr[k]));
      }

      // sort the entries by column
      std::sort(entries.begin(), entries.end(), [](const std::tuple<int, int, double> &a, const std::tuple<int, int, double> &b) {
        return std::get<1>(a) < std::get<1>(b);
      });

      derived_vcsc_mat->mat = new IVSparse::VCSC<double, int>(entries, nrow, ncol, nnz);
      vcsc_mat = derived_vcsc_mat;

      return;
    }

    if (columnMajor) {
      DerivedVCSC<double, int> *derived_vcsc_mat = new DerivedVCSC<double, int>();
      derived_vcsc_mat->mat = new IVSparse::VCSC<double, int>(x_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
      vcsc_mat = derived_vcsc_mat;
    } else {
      DerivedVCSC<double, int, 0> *derived_vcsc_mat = new DerivedVCSC<double, int, 0>();
      derived_vcsc_mat->mat = new IVSparse::VCSC<double, int, 0>(x_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
      vcsc_mat = derived_vcsc_mat;
    }
  }

  // Constructor to convert a dgCMatrix to VCSC with specified data and index types
  VCSC(const S4 &mat, std::string data_type, std::string index_type) {

    // create the vectors
    IntegerVector i, p, Dim;
    NumericVector x;
    List Dimnames;
    bool columnMajor = true;

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

    if (data_type == "int" && index_type == "int") {
      DerivedVCSC<int, int> *derived_vcsc_mat = new DerivedVCSC<int, int>();
      std::vector<int> x2(x.begin(), x.end());     // values
      int *x2_ptr = &x2[0];
      derived_vcsc_mat->mat = new IVSparse::VCSC<int, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
      vcsc_mat = derived_vcsc_mat;
    } 

    std::vector<uint64_t> double_indices;
    std::vector<uint64_t> double_pointers;
    uint64_t *i2_ptr;
    uint64_t *p2_ptr;
    int index_byte_size = 4;
    if (index_type == "long" || index_type == "uint64_t") {
      index_byte_size = 8;
      double_indices = std::vector<uint64_t>(i.begin(), i.end());
      double_pointers = std::vector<uint64_t>(p.begin(), p.end());
      i2_ptr = &double_indices[0];
      p2_ptr = &double_pointers[0];
    }

    if (data_type == "int") { // ----------------- int ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<int, int> *derived_vcsc_mat = new DerivedVCSC<int, int>();
        std::vector<int> x2(x.begin(), x.end());     // values
        int *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<int, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<int, uint64_t> *derived_vcsc_mat = new DerivedVCSC<int, uint64_t>();
        std::vector<int> x2(x.begin(), x.end());     // values
        int *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<int, uint64_t>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else if (data_type == "double") { // ----------------- double ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<double, int> *derived_vcsc_mat = new DerivedVCSC<double, int>();
        derived_vcsc_mat->mat = new IVSparse::VCSC<double, int>(x_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<double, uint64_t> *derived_vcsc_mat = new DerivedVCSC<double, uint64_t>();
        derived_vcsc_mat->mat = new IVSparse::VCSC<double, uint64_t>(x_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else if (data_type == "float") { // ----------------- float ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<float, int> *derived_vcsc_mat = new DerivedVCSC<float, int>();
        std::vector<float> x2(x.begin(), x.end());     // values
        float *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<float, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<float, uint64_t> *derived_vcsc_mat = new DerivedVCSC<float, uint64_t>();
        std::vector<float> x2(x.begin(), x.end());     // values
        float *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<float, uint64_t>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else if (data_type == "long") { //----------------- long ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<long, int> *derived_vcsc_mat = new DerivedVCSC<long, int>();
        std::vector<long> x2(x.begin(), x.end());     // values
        long *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<long, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<long, uint64_t> *derived_vcsc_mat = new DerivedVCSC<long, uint64_t>();
        std::vector<long> x2(x.begin(), x.end());     // values
        long *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<long, uint64_t>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else if (data_type == "short") { //----------------- short ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<short, int> *derived_vcsc_mat = new DerivedVCSC<short, int>();
        std::vector<short> x2(x.begin(), x.end());     // values
        short *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<short, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<short, uint64_t> *derived_vcsc_mat = new DerivedVCSC<short, uint64_t>();
        std::vector<short> x2(x.begin(), x.end());     // values
        short *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<short, uint64_t>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else if (data_type == "uint8_t") { //----------------- uint8_t ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<uint8_t, int> *derived_vcsc_mat = new DerivedVCSC<uint8_t, int>();
        std::vector<uint8_t> x2(x.begin(), x.end());     // values
        uint8_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint8_t, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<uint8_t, uint64_t> *derived_vcsc_mat = new DerivedVCSC<uint8_t, uint64_t>();
        std::vector<uint8_t> x2(x.begin(), x.end());     // values
        uint8_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint8_t, uint64_t>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else if (data_type == "uint16_t") { //----------------- uint16_t ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<uint16_t, int> *derived_vcsc_mat = new DerivedVCSC<uint16_t, int>();
        std::vector<uint16_t> x2(x.begin(), x.end());     // values
        uint16_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint16_t, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<uint16_t, uint64_t> *derived_vcsc_mat = new DerivedVCSC<uint16_t, uint64_t>();
        std::vector<uint16_t> x2(x.begin(), x.end());     // values
        uint16_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint16_t, uint64_t>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else if (data_type == "uint16_t") { //----------------- uint16_t ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<uint16_t, int> *derived_vcsc_mat = new DerivedVCSC<uint16_t, int>();
        std::vector<uint16_t> x2(x.begin(), x.end());     // values
        uint16_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint16_t, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<uint16_t, uint64_t> *derived_vcsc_mat = new DerivedVCSC<uint16_t, uint64_t>();
        std::vector<uint16_t> x2(x.begin(), x.end());     // values
        uint16_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint16_t, uint64_t>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else if (data_type == "uint32_t") { //----------------- uint32_t ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<uint32_t, int> *derived_vcsc_mat = new DerivedVCSC<uint32_t, int>();
        std::vector<uint32_t> x2(x.begin(), x.end());     // values
        uint32_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint32_t, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<uint32_t, uint64_t> *derived_vcsc_mat = new DerivedVCSC<uint32_t, uint64_t>();
        std::vector<uint32_t> x2(x.begin(), x.end());     // values
        uint32_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint32_t, uint64_t>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else if (data_type == "uint64_t") { //----------------- uint64_t ----------------- //
      
      if (index_byte_size == 4) {
        DerivedVCSC<uint64_t, int> *derived_vcsc_mat = new DerivedVCSC<uint64_t, int>();
        std::vector<uint64_t> x2(x.begin(), x.end());     // values
        uint64_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint64_t, int>(x2_ptr, i_ptr, p_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      } else {
        DerivedVCSC<uint64_t, uint64_t> *derived_vcsc_mat = new DerivedVCSC<uint64_t, uint64_t>();
        std::vector<uint64_t> x2(x.begin(), x.end());     // values
        uint64_t *x2_ptr = &x2[0];
        derived_vcsc_mat->mat = new IVSparse::VCSC<uint64_t, uint64_t>(x2_ptr, i2_ptr, p2_ptr, nrow, ncol, nnz);
        vcsc_mat = derived_vcsc_mat;
      }

    } else {
      
      stop("Invalid data type");

    }
  }

  // Constructor for file input
  // VCSC(std::string filename) {}

  // Copy constructor
  // VCSC(const VCSC &vcsc) {
  //   vcsc_mat = vcsc.vcsc_mat;
  // }

  // IVCSC constructor
  // VCSC(const IVCSC &ivcsc) {}

  //-------------------------------------------------//

  // ------- Destructor ------- //
  ~VCSC() { delete vcsc_mat; }

  /****************************************************
   *                                                   *
   *                                                   *
   *                   Converters                      *
   *                                                   *
   *                                                   *
   *****************************************************/

  // to dgCMatrix
  // S4 to_dgCMatrix() {}

  // // to dgRMatrix
  // S4 to_dgRMatrix() {}

  // // to dgTMatrix
  // S4 to_dgTMatrix() {}

  // to ivcsc
  // IVCSC to_ivcsc() {}

  /****************************************************
   *                                                   *
   *                                                   *
   *                      Getters                      *
   *                                                   *
   *                                                   *
   *****************************************************/

  // get the value at position (i, j)
  double coeff(int i, int j) { return vcsc_mat->coeff(i, j); }

  // rows
  int rows() { return vcsc_mat->rows(); }

  // cols
  int cols() { return vcsc_mat->cols(); }

  // nnz
  int nnz() { return vcsc_mat->nnz(); }

  // innerdim
  int innerdim() { return vcsc_mat->innerdim(); }

  // outerdim
  int outerdim() { return vcsc_mat->outerdim(); }

  // bytesize
  int bytesize() { return vcsc_mat->bytesize(); }

  /****************************************************
   *                                                   *
   *                                                   *
   *               Matrix Manipulation                 *
   *                                                   *
   *                                                   *
   *****************************************************/

  // in place transpose
  void transpose() { vcsc_mat->transpose();}

  // // slice
  // VCSC slice(int start, int end) {}

  // append
  void append(VCSC &other) { vcsc_mat->append(other); } 

  // matrix multiplication
  NumericMatrix mult(NumericMatrix dense_mat) { return vcsc_mat->mult(dense_mat); }

  /*****************************************************
   *                                                   *
   *                                                   *
   *                     Operators                     *
   *                                                   *
   *                                                   *
   *****************************************************/

  // scale
  // VCSC operator*(double scalar) {}

  // // spmm
  // VCSC operator*(NumericMatrix dense_mat) {}

  // // equality
  // bool operator==(VCSC vcsc) {}

  // // inequality
  // bool operator!=(VCSC vcsc) {}

  // // assignment
  // VCSC operator=(VCSC vcsc) {}

  /*****************************************************
   *                                                   *
   *                                                   *
   *                       Misc                        *
   *                                                   *
   *                                                   *
   *****************************************************/

  // print
  // void print() {}

  // // save
  // void save(std::string filename) {}

  // // load
  // VCSC load(std::string filename) {}

  //------------------------------------------------------

//  private:
  BaseVCSC *vcsc_mat;
};



// // DerivedVCSC methods that have circular dependencies
template <typename T, typename U, bool columnMajor>
void DerivedVCSC<T, U, columnMajor>::append(VCSC& other) {
  IVSparse::VCSC<T, U, columnMajor> *other_mat;
  // other_mat = other.vcsc_mat->mat;
  
  // get number of columns of other
  int other_cols = other.cols();
  
  //print number of added columns
  Rprintf("Adding %d columns\n", other_cols);
}




// Expose the VCSC class
RCPP_MODULE(vcsc) {
  class_<VCSC>("VCSC")
    .constructor()
    .constructor<S4>()
    .constructor<S4, std::string, std::string>()
    // .constructor<std::string>()
    // .constructor<VCSC>()
    // // .constructor<IVCSC>()
    // // Add methods here
    // // converter methods
    // .method("to_dgCMatrix", &VCSC::to_dgCMatrix)
    // .method("to_dgRMatrix", &VCSC::to_dgRMatrix)
    // .method("to_dgTMatrix", &VCSC::to_dgTMatrix)
    // // getters
    .method("coeff", &VCSC::coeff)
    .method("rows", &VCSC::rows)
    .method("cols", &VCSC::cols)
    .method("nnz", &VCSC::nnz)
    .method("innerdim", &VCSC::innerdim)
    .method("outerdim", &VCSC::outerdim)
    .method("bytesize", &VCSC::bytesize)
    // // matrix manipulation
    .method("transpose", &VCSC::transpose)
    // .method("slice", &VCSC::slice)
    .method("append", &VCSC::append)
    .method("mult", &VCSC::mult)
    // // operators
    // .method("operator*", &VCSC::operator*)
    // .method("operator*", &VCSC::operator*)
    // .method("operator==", &VCSC::operator==)
    // .method("operator!=", &VCSC::operator!=)
    // .method("operator=", &VCSC::operator=)
    // // misc
    // .method("print", &VCSC::print)
    // .method("save", &VCSC::save)
    // .method("load", &VCSC::load)
    ;
}
