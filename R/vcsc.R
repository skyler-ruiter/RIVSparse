library(R6)

VCSC <- R6::R6Class("VCSC",
  public = list(
    vcsc_instance = NULL,
    value_type = NULL,
    index_type = NULL,
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
        ptr = self$vcsc_instance$get_vcsc_xptr() # get ivsparse ptr
        new_vcsc_instance$set_vcsc_ptr(ptr) # make new ivsparse object and set ptr
        temp = VCSC$new() # make new VCSC object
        temp$vcsc_instance <- new_vcsc_instance # set vcsc_instance to new_vcsc_instance
        return(temp)
      } else {
        return(VCSC$new(self))
      }
    },
    coeff = function(i, j) {
      self$vcsc_instance$coeff(i, j)
    },
    scale = function(alpha) {
      self$vcsc_instance$scale(alpha)
    },
    print = function() {
      self$vcsc_instance$print()
    },
    append = function(other) {
      self$vcsc_instance$append(other$vcsc_instance$get_vcsc_xptr())
    },
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
    test2 = function(mat) {
      print(mat$vcsc_instance)
      self$vcsc_instance$test2(mat)
    },
    test_return = function(mat) {
      return(self$vcsc_instance$test_return(mat))
    },
    get_vcsc_xptr = function() {
      return(self$vcsc_instance$get_vcsc_xptr())
    },
    set_vcsc_ptr = function(mat) {
      ptr = mat$get_vcsc_xptr()
      self$vcsc_instance$set_vcsc_ptr(ptr)
    }
  )
)

# setClass("VCSC",
#   representation = representation(
#     ref = "VCSC_ref"
#   )
# )