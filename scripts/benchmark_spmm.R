library(RIVSparse)
library(Rcpp)
library(Matrix)
library(microbenchmark)

n <- 5000

dense <- matrix(rnorm(n*n), ncol=n) # dense matrix of doubles

sparse <- rsparsematrix(n, n, 0.1) # random (mostly unique) matrix of doubles

vcsc1 <- VCSC$new(sparse) # use doubles and int indices

orig_spmm <- function(sparse, dense) {sparse %*% dense}

new_spmm <- function(sparse, dense) {sparse$mult(dense)}

bench <- microbenchmark(
  orig_spmm(sparse, dense),
  new_spmm(vcsc1, dense),
  times = 3
)

print(bench)
