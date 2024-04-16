library(R6)

VCSC <- R6::R6Class("VCSC",
  public = list(
    vcsc_instance = NULL,
    initialize = function(A, val_t = "double", idx_t = "int") {
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
    slice = function(i, j) {
      # get and return the slice
      self$vcsc_instance$slice(i, j)
    }
  )
)

# setClass("VCSC",
#   representation = representation(
#     ref = "VCSC_ref"
#   )
# )