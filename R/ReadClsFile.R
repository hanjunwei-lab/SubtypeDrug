##' ReadClsFile
##'
##'
##' @title SubtypeDrug internal function
##' @description These are function read sample label file (.cls format).
##' @param file Input sample subtype class vector file in CLS format.
##' @return a list
##' @examples
##' Subtype<-system.file("extdata", "Subtype_labels.cls", package = "SubtypeDrug")
##' x<-ReadClsFile(Subtype)
##' @author Xudong Han,
##' Junwei Han,
##' Chonghui Liu
##' @export
ReadClsFile<-function(file) {
  cls.cont <- readLines(file)
  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  t <- rev(table(class.list))
  l <- length(t)
  phen <- vector(length=l, mode="character")
  for (i in 1:l) {
    phen[i] <- noquote(names(t)[i])
  }
  return(list(phen = phen, class.labes = class.list))
}
