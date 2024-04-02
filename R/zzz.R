.onAttach <- function(libname, pkgname) {

  # ensure rng seed is initialized
  if (!exists(".Random.seed")) set.seed(NULL)

  # set threads
  if (is.null(getOption("RIVSparse.threads"))) {
    options(RIVSparse.threads = 0)
  }

  if (is.null(getOption("RIVSparse.verbose"))) {
    options(RIVSparse.verbose = FALSE)
  }

  if (is.null(getOption("RIVSparse.debug"))) {
    options(RIVSparse.debug = TRUE)
  }

  msg <- "RIVSparse v0.1.0 using 'options(RIVSparse.threads = "
  threads <- getOption("RIVSparse.threads")
  if (threads == 0) {
    msg <- c(msg, "0)' (all available threads)")
  } else {
    msg <- c(msg, threads, ")")
  }
  msg <- c(msg, ", 'options(RIVSparse.verbose = ")
  if (getOption("RIVSparse.verbose")) {
    msg <- c(msg, "TRUE)'\n")
  } else {
    msg <- c(msg, "FALSE)'")
  }

  packageStartupMessage(msg)
  invisible()

}