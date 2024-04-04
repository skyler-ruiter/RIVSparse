library(RIVSparse)
library(Rcpp)
library(Matrix)
library(microbenchmark)

n <- 5000
low <- 10
high <- 20


dense <- matrix(rnorm(n * n), ncol = n) # dense matrix of doubles

sparse <- rsparsematrix(n, n, 0.1) # random (mostly unique) matrix of doubles
vcsc1 <- VCSC$new(sparse) # use doubles and int indices

sparse_redundant <- rsparsematrix(n, n, 0.1, rand.x = function(n) rpois(2, 1) + 1)
vcsc2 <- VCSC$new(sparse_redundant, "int", "int")

orig_spmm <- function(sparse, dense) {
  sparse %*% dense
}

new_spmm <- function(sparse, dense) {
  sparse$mult(dense)
}

bench <- microbenchmark(
  orig_spmm(sparse, dense),
  new_spmm(vcsc1, dense),
  times = 3
)

bench2 <- microbenchmark(
  orig_spmm(sparse_redundant, dense),
  new_spmm(vcsc2, dense),
  times = 3
)

print("Unique Sparse")
print(bench)

print("50% redundant Sparse")
print(bench2)
