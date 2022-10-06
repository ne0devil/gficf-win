.onAttach <- function(libname, pkgname) {
  require(Matrix)
  packageStartupMessage(packageVersion("gficf"))
  invisible()
}