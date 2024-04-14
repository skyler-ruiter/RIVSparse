ignore_me <- setMethod("show", "Rcpp_VCSC", function(object) {
  cat("\nMatrix of type RIVSparse::VCSC\n")

  for (i in 0:(object$outerdim() - 1)) {
    cat("\n")
    for (j in 0:(object$innerdim() - 1)) {
      cat(object$coeff(i, j), " ")
    }
  }
  cat("\n")

})