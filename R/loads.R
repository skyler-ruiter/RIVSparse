library(Rcpp)
library(methods)

# load the module into the namespace for end user
Rcpp::loadModule("yada", TRUE)

# load vcsc module
Rcpp::loadModule("vcsc", TRUE)
