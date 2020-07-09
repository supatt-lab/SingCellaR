python.fa2 = NULL
python.nx  = NULL
python.np  = NULL

.onLoad = function(libname, pkgname) {
  if (suppressWarnings(suppressMessages(requireNamespace("reticulate")))) {
    has.pkg.fa2 = reticulate::py_module_available("fa2")
    has.pkg.nx = reticulate::py_module_available("networkx")
    has.pkg.np = reticulate::py_module_available("numpy")
    
    if (has.pkg.fa2) {
      python.fa2 <<- reticulate::import("fa2", delay_load=TRUE)
    }
    if (has.pkg.nx) {
      python.nx <<- reticulate::import("networkx", delay_load=TRUE)
    }
    if (has.pkg.np) {
      python.np <<- reticulate::import("numpy", convert=FALSE,delay_load=TRUE)
    }
  }
}