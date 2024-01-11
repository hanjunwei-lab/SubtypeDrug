##' isPackageLoaded
##'
##'
##' @title SubtypeDrug internal function
##' @description Determine if the package is loaded. If the package is not loaded,
##' the program will prompt the user.
##' @param name A string. The name of the R package which determines whether it is loaded.
##' @return A string, TRUE or FALSE.
##' @author Xudong Han,
##' Junwei Han,
##' Chonghui Liu
##' @examples
##' isPackageLoaded("pheatmap")
##' @export
isPackageLoaded <- function(name) {
  ## Purpose: is package 'name' loaded?
  ## --------------------------------------------------
  (paste("package:", name, sep="") %in% search()) ||
    (name %in% loadedNamespaces())
}
