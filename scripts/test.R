library(RIVSparse)
library(Matrix)

# make a random sparse 10x10 matrix of csc, csr, and coo
set.seed(1)
n <- 10
A <- rsparsematrix(n, n, 1, repr="C") # csc
B <- rsparsematrix(n, n, 0.5, repr="R") # csr
C <- rsparsematrix(n, n, 0.5, repr="T") # coo

# vcsc_A <- VCSC$new(A)
# vcsc_B <- VCSC$new(B)
# vcsc_C <- VCSC$new(C)

A2 <- rsparsematrix(n, n, 0.5, rand.x=function(n) rpois(n, 1) + 1)

# temp <- new(VCSC_INT_INT, A2)

# make an R object of VCSC 
vcsc_A <- VCSC$new(A2, "int", "int")
vcsc2_A <- VCSC$new(A2, "int", "int")

vcsc_A

