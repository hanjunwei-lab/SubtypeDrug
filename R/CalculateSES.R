##' CalculateSES
##'
##'
##' @title SubtypeDrug internal function
##' @description Calculate subpathway enrichment score.
##' @param labels.list A vector of 0 and 1.
##' @param correl.vector A vector. The weight value used to calculate the enrichment score.
##' @return A vector.
##' @author Xudong Han,
##' Junwei Han,
##' Chonghui Liu
##' @examples
##' x<-CalculateSES(sample(c(0,1),10,replace = TRUE),c(1:10))
##' @export
CalculateSES<-function(labels.list,correl.vector = NULL){
  tag.indicator <- labels.list
  no.tag.indicator <- 1 - tag.indicator
  N <- length(labels.list)
  Nh <- length(labels.list[labels.list==1])
  Nm <-  N - Nh
  correl.vector <- abs(correl.vector)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)

  max.ES <- max(RES)
  min.ES <- min(RES)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  return(ES)
}

