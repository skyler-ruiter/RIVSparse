library(R6)

VCSC <- R6Class(
  classname = "VCSC",
  private = list(
    cpp_ptr = NULL
  ),
  initialize = function(data) {
    self$cpp_ptr <- create_vcsc(data)
  }
)

