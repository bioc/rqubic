.onLoad <- function(lib, pkg) {
  ##  library.dynam("rqubic", pkg, lib)
  .Call("RQUBIC_init", pkg)
}
