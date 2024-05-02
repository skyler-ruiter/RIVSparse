library(RIVSparse)
library(Matrix)

# make a random sparse 10x10 matrix of csc, csr, and coo

set.seed(1)
n <- 10

A <- rsparsematrix(n, n, 1, repr="C") # csc
B <- rsparsematrix(n, n, 0.5, repr="R") # csr
C <- rsparsematrix(n, n, 0.5, repr="T") # coo

# A
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

# A2

temp <- new(VCSC_INT_INT, A2)


# make an R object of VCSC 
vcsc_A <- VCSC$new(A2, "int", "int")
vcsc2_A <- VCSC$new(A2, "int", "int")

# vcsc_B <- VCSC$new(A)

# run coeff on (1, 1)
# vcsc_A$coeff(1, 1)
# vcsc_B$coeff(1, 1)

# append vcsc_A to vcsc2_A
# vcsc2_A$append(vcsc_A)

# print("hi")

# # scale vcsc_A by 2
# vcsc_A$scale(2)

# # print vcsc_A
# vcsc_A$print()

vcsc_A$test2(vcsc2_A)

# test <- VCSC$new(A2)

# res <- vcsc_A$test_return(vcsc2_A)

# res

# vcsc_A$test_ptr(vcsc2_A)

# print the index and value types of vcsc_A
print(vcsc_A$value_type)
print(vcsc_A$index_type)

print(vcsc_A$get_vcsc_xptr())

# vcsc2_A$print()

# copy_vcsc <- vcsc_A$clone()
# deep_copy_vcsc <- vcsc_A$clone(deep = TRUE)


# vcsc2_A$print()

# vcsc2_A$set_vcsc_ptr(vcsc_A)

# vcsc2_A$print()


# test copy
vcsc_new <- vcsc_A$copy(deep = TRUE)

vcsc_A$vcsc_instance

vcsc_new$vcsc_instance

vcsc_A$scale(2)

vcsc_A$print()

vcsc_new$print()

vcsc_new$append(vcsc_A)

vcsc_new$print()

vcsc_slice <- vcsc_A$slice(1, 5)

vcsc_slice$print()