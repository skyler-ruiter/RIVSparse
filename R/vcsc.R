library(R6)

VCSC <- R6::R6Class("VCSC",
  public = list(
    vcsc_instance = NULL,
    value_type = NULL,
    index_type = NULL,
    initialize = function(A, val_t = "double", idx_t = "int") {

      self$value_type <- val_t
      self$index_type <- idx_t

      # if A is a VCSC object
      if (inherits(A, "VCSC")) {
        # check data types match
        if (A$value_type != val_t || A$index_type != idx_t) {
          stop("Data types do not match")
        }
        self$vcsc_instance <- A$vcsc_instance
        printf("Hi from VCSC copy constructor section\n")
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
    coeff = function(i, j) {
      self$vcsc_instance$coeff(i, j)
    },
    append = function(other) {
      self$vcsc_instance$append(other)
    },
    scale = function(alpha) {
      self$vcsc_instance$scale(alpha)
    },
    print = function() {
      self$vcsc_instance$print()
    },
    test2 = function(mat) {
      print(mat$vcsc_instance)
      self$vcsc_instance$test2(mat)
    }
  )
)

# setClass("VCSC",
#   representation = representation(
#     ref = "VCSC_ref"
#   )
# )