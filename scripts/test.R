library(RIVSparse)
library(Matrix)

# make a random sparse 10x10 matrix of csc, csr, and coo

set.seed(1)
n <- 10

A <- rsparsematrix(n, n, 1, repr="C") # csc
B <- rsparsematrix(n, n, 0.5, repr="R") # csr
C <- rsparsematrix(n, n, 0.5, repr="T") # coo

A
# B
# C

# vcsc_A <- VCSC$new(A)
# vcsc_B <- VCSC$new(B)
# vcsc_C <- VCSC$new(C)

# transpose vcsc_A
# vcsc_A$transpose()

# str(A)
# str(B)
# str(C)

A2 <- rsparsematrix(n, n, 0.5, rand.x=function(n) rpois(n, 1) + 1)

A2

temp <- new(VCSC_INT_INT, A2)

print("hi")

# make an R object of VCSC 
vcsc_A <- VCSC$new(A2, "int", "int")
vcsc2_A <- VCSC$new(A2, "int", "int")

vcsc_B <- VCSC$new(A)

# run coeff on (1, 1)
vcsc_A$coeff(1, 1)
vcsc_B$coeff(1, 1)

print("hi")

isS4(vcsc_A)

print(vcsc_A)

# append vcsc_A to vcsc2_A
# vcsc2_A$append(vcsc_A)

# print("hi")

# scale vcsc_A by 2
vcsc_A$scale(2)

# print vcsc_A
vcsc_A$print()