##' AccumulateNormal
##'
##'
##' @title SubtypeDrug internal function
##' @description Infering patient-specific subpathway activity profiles.
##' @param x_matrix A subpathway activity profile. rows are subpathwyas, columns are samples.
##' @param control_index A vector. In the sample of the subpathway activity profile, the position of control samples.
##' @return A matrix.
##' @importFrom stats var
##' @importFrom stats sd
##' @author Xudong Han,
##' Junwei Han,
##' Chonghui Liu
##' @examples
##' x<-matrix(c(1:10),ncol = 5)
##' x1<-AccumulateNormal(x,c(3,5))
##' @export
AccumulateNormal<-function(x_matrix,control_index){
  control_matrix<-x_matrix[,control_index]
  control_mean<-apply(control_matrix, 1, mean)
  control_sd<-apply(control_matrix, 1, sd)
  spw_matrix<-x_matrix[,-control_index]
  spw_matrix<-apply(spw_matrix, 2, function(x){
    x<-(x-control_mean)/control_sd
    return(x)
  })
  return(spw_matrix)
}
