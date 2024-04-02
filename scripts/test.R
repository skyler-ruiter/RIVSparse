library(RIVSparse)
library(Matrix)

n <- 10
a <- rsparsematrix(n, n, 0.5, rand.x=function(n) rpois(n, 1) + 1)

a

mat <- VCSC$new(a)


mat$coeff(1,1)

dense <- matrix(sample.int(50), 10, 10)

mat$mult(dense)

a %*% dense

