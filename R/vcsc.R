library(R6)

VCSC <- R6::R6Class("VCSC",
  public = list(

    ###* public fields *###
    vcsc_instance = NULL,
    value_type = NULL,
    index_type = NULL,

    ###* Constructors and Copy *###
    initialize = function(A = NULL, val_t = "double", idx_t = "int") {

      self$value_type <- val_t
      self$index_type <- idx_t

      if (is.null(A)) {
        # empty initialization
        return()
      }

      # if A is a VCSC object
      if (inherits(A, "VCSC")) {
        print("Warning: Shallow copy of VCSC object")
        self$vcsc_instance <- A$vcsc_instance
        return()
      }

      if (val_t == "int" && idx_t == "int") {
        self$vcsc_instance <- new(VCSC_INT_INT, A)
      } else if (val_t == "double" && idx_t == "int") {
        self$vcsc_instance <- new(VCSC_DOUBLE_INT, A)
      } else {
        stop("Invalid type")
      }

    },

    # clone method for deep and shallow copy
    copy = function(deep = FALSE) {
      val_t <- self$value_type
      idx_t <- self$index_type
      new_vcsc_instance <- NULL
      if (deep) {
        # resolve template types
        if (val_t == "int" && idx_t == "int") {
          # make a deep copy of VCSC_INT_INT
          new_vcsc_instance <- new(VCSC_INT_INT)
        } else if (val_t == "double" && idx_t == "int") {
          # make a deep copy of VCSC_DOUBLE_INT
          new_vcsc_instance <- new(VCSC_DOUBLE_INT)
        }
        ptr <- self$vcsc_instance$get_vcsc_xptr() # get ivsparse ptr
        # make new ivsparse object and set ptr
        new_vcsc_instance$set_vcsc_ptr(ptr)
        temp <- VCSC$new() # make new VCSC object
        # set vcsc_instance to new_vcsc_instance
        temp$vcsc_instance <- new_vcsc_instance
        return(temp)
      } else {
        return(VCSC$new(self))
      }
    },

    ###* Getters *###

    # get value at (i, j)
    coeff = function(i, j) {
      self$vcsc_instance$coeff(i, j)
    },

    # get number of rows
    rows = function() {
      return(self$vcsc_instance$rows())
    },

    # get number of columns
    cols = function() {
      return(self$vcsc_instance$cols())
    },

    # get inner dimension
    innerDim = function() {
      return(self$vcsc_instance$innerDim())
    },

    # get outer dimension
    outerDim = function() {
      return(self$vcsc_instance$outerDim())
    },

    # get nnz
    nnz = function() {
      return(self$vcsc_instance$nnz())
    },

    # get byte size
    byteSize = function() {
      return(self$vcsc_instance$byteSize())
    },

    # get major order
    getColumnMajor = function() {
      return(self$vcsc_instance$getColumnMajor())
    },

    # get the unique values for the column
    getUniqueValues = function(j) {
      return(self$vcsc_instance$getValues(j))
    },

    # get the counts for the column
    getCounts = function(j) {
      return(self$vcsc_instance$getCounts(j))
    },

    # get the indices for the column
    getIndices = function(j) {
      return(self$vcsc_instance$getIndices(j))
    },

    # get the number of unique values for the column
    getNumUniqueVals = function(j) {
      return(self$vcsc_instance$getNumUniqueVals(j))
    },

    # get num indices for the column
    getNumIndices = function(j) {
      return(self$vcsc_instance$getNumIndices(j))
    },

    ###* Converters *###

    ###* Calculations *###

    ###* Utilities *###

    # print method
    print = function() {
      self$vcsc_instance$print()
    },

    ###* Matrix Manipulation *###

    # append other VCSC object to self
    append = function(other) {
      self$vcsc_instance$append(other$vcsc_instance$get_vcsc_xptr())
    },

    # slice VCSC object
    slice = function(start, end) {
      ptr = self$vcsc_instance$slice(start, end)
      val_t <- self$value_type
      idx_t <- self$index_type
      new_vcsc_instance <- NULL
      if (val_t == "int" && idx_t == "int") {
        new_vcsc_instance <- new(VCSC_INT_INT)
      } else if (val_t == "double" && idx_t == "int") {
        new_vcsc_instance <- new(VCSC_DOUBLE_INT)
      }
      new_vcsc_instance$set_vcsc_ptr(ptr)
      temp = VCSC$new()
      temp$vcsc_instance <- new_vcsc_instance
      return(temp)
    },

    ###* Matrix Operations *###

    # in place scalar multiplication
    scale = function(alpha) {
      self$vcsc_instance$scale(alpha)
    },

    ###* Operator Overloads *###

    ###* R6 Methods *###

    # get underlying ivsparse xptr
    get_vcsc_xptr = function() {
      return(self$vcsc_instance$get_vcsc_xptr())
    },

    # set underlying ivsparse xptr
    set_vcsc_ptr = function(mat) {
      ptr = mat$get_vcsc_xptr()
      self$vcsc_instance$set_vcsc_ptr(ptr)
    }

    ###* Testing *###

    # test2 = function(mat) {
    #   print(mat$vcsc_instance)
    #   self$vcsc_instance$test2(mat)
    # },
    # test_return = function(mat) {
    #   return(self$vcsc_instance$test_return(mat))
    # },

  )
)
