.onAttach <- function(libname, pkgname) {
  require(Matrix)
  packageStartupMessage("gficf v0.5.0")
}