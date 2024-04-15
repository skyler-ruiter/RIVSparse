library(RIVSparse)
library(Matrix)

# make a random sparse 10x10 matrix of csc, csr, and coo

set.seed(1)
n <- 10

A <- rsparsematrix(n, n, 0.5, repr="C") # csc
B <- rsparsematrix(n, n, 0.5, repr="R") # csr
C <- rsparsematrix(n, n, 0.5, repr="T") # coo

# A
# B
# C

vcsc_A <- VCSC$new(A)
vcsc_B <- VCSC$new(B)
vcsc_C <- VCSC$new(C)

# transpose vcsc_A
vcsc_A$transpose()

# str(A)
# str(B)
# str(C)


A2 <- rsparsematrix(n, n, 0.5, rand.x=function(n) rpois(n, 1) + 1)

# test coeff, rows, cols, nnz, inner, outer, and bytesize
# vcsc_A$coeff(1,1)
# vcsc_A$coeff(1,2)
# vcsc_A$rows()
# vcsc_A$cols()
# vcsc_A$nnz()
# vcsc_A$innerdim()
# vcsc_A$outerdim()
# vcsc_A$bytesize()

vcsc_A2 <- VCSC$new(A2, "int", "int")
vcsc_A2_copy <- new(VCSC, A2, "int", "int")

vcsc_A2

# make a reference to vcsc_A2
vcsc_A2_ref <- vcsc_A2

# test append
# vcsc_A2$append(vcsc_A2_copy)


vcsc_A2

# make a 1000x1000 matrix
# n <- 1000
# A2 <- rsparsematrix(n, n, 0.5, rand.x=function(n) rpois(n, 1) + 1)
# temp <- VCSC$new(A2, "int", "int")
# temp

print("Testing .print()")
vcsc_A2$print()


# A
# A2

# vcsc_int <- VCSC$new(A2, "int", "int")
# vcsc_int$coeff(1,1)

# vcsc_double <- VCSC$new(A, "double", "int")
# vcsc_double$coeff(1,1)

# vcsc_float <- VCSC$new(A, "float", "int")
# vcsc_float$coeff(1,1)

# vcsc_long <- VCSC$new(A2, "long", "int")
# vcsc_long$coeff(1,1)

# vcsc_short <- VCSC$new(A2, "short", "int")
# vcsc_short$coeff(1,1)

# vcsc_uint8 <- VCSC$new(A2, "uint8_t", "int")
# vcsc_uint8$coeff(1,1)

# vcsc_uint16 <- VCSC$new(A2, "uint16_t", "int")
# vcsc_uint16$coeff(1,1)

# vcsc_uint32 <- VCSC$new(A2, "uint32_t", "int")
# vcsc_uint32$coeff(1,1)

# vcsc_uint64 <- VCSC$new(A2, "uint64_t", "int")
# vcsc_uint64$coeff(1,1)

# vcsc_int_uint64 <- VCSC$new(A2, "int", "uint64_t")
# vcsc_int_uint64$coeff(1,1)