ignore_me <- setMethod("show", "Rcpp_VCSC", function(object) {
  cat("\nMatrix of type RIVSparse::VCSC\n")

  # if more than 100 rows/cols don't print the whole matrix
  # just print a summary with the dimensions, nnz and the first 10 rows/cols
  if (object$outerdim() > 100 || object$innerdim() > 100) {
    cat("Dimensions: ", object$outerdim(), "x", object$innerdim(), "\n")
    cat("Non-zero elements: ", object$nnz(), "\n")
    cat("First 10 rows/cols:\n")
    for (i in 0:9) {
      cat("\n")
      for (j in 0:9) {
        cat(object$coeff(i, j), " ")
      }
    }
    cat("\n")
    return()
  } else {
    for (i in 0:(object$outerdim() - 1)) {
      cat("\n")
      for (j in 0:(object$innerdim() - 1)) {
        cat(object$coeff(i, j), " ")
      }
    }
    cat("\n")
  }

})