##' getDrugSpw
##'
##'
##' @title SubtypeDrug internal function
##' @description According to the parameters set by the user, the up-regulatory
##' and down-regulatory subpathway data of drug is obtained.
##' @param drug_target_data A list. A list stores a collection of drug up- and down-regulated subpathways.
##' @param spw_matrix_rnames A vector. A vector consisting of row names of subpathway activity profile.
##' @param drug.P.value.threshold A value. According to the threshold of the significant P value set by parameter `drug.p.val.threshold`, the drug up-regulation and down-regulatory subpathways were screened.
##' @param drug.min.sz  A numeric. The drug regulated subpathways intersects with the subpathways in the subpathway activity profile. Then drugs with less than `drug.spw.min.sz` up- or down-regulated subpathways are removed.
##' @param drug.max.sz A numeric. Similar to parameter `drug.spw.min.sz`, drugs with more than `drug.spw.max.sz` up- or down-regulated subpathways are removed.
##' @return a list.
##' @examples
##' require(GSVA)
##' Geneexp<-get("Geneexp")
##' UserGS<-get("UserGS")
##' UserDS<-get("UserDS")
##' spw_matrix<-gsva(Geneexp,UserGS,verbose=FALSE)
##' x<-getDrugSpw(UserDS,row.names(spw_matrix),0.05,1,100)
##' @author Xudong Han,
##' Junwei Han,
##' Chonghui Liu
##' @export
getDrugSpw<-function(drug_target_data,spw_matrix_rnames,drug.P.value.threshold,drug.min.sz,drug.max.sz){
  isdf<-is.data.frame(drug_target_data[[1]])
  if(isdf==TRUE){
   Target_upregulation<-lapply(drug_target_data, function(x){
              screen_index<-which(x$ES>0&x$P_value<=drug.P.value.threshold)
              drug_up<-x$SubpathwayId[screen_index]
              return(drug_up)
              })
   Target_downregulation<-lapply(drug_target_data, function(x){
              screen_index<-which(x$ES<0&x$P_value<=drug.P.value.threshold)
              drug_down<-x$SubpathwayId[screen_index]
              return(drug_down)
              })
  }else{
    Target_upregulation<-list()
    Target_downregulation<-list()
    for(i in 1:length(drug_target_data)){
      Target_upregulation[[i]]<-drug_target_data[[i]]$Target_upregulation
      Target_downregulation[[i]]<-drug_target_data[[i]]$Target_downregulation
    }
    names(Target_upregulation)<-names(drug_target_data)
    names(Target_downregulation)<-names(drug_target_data)
  }
   Target_upregulation<- lapply(Target_upregulation,function(x ,y){
              x1<-x[which(is.element(x,y)==TRUE)]
              return(x1)
              },spw_matrix_rnames)
   Target_downregulation<- lapply(Target_downregulation,function(x ,y){
              x1<-x[which(is.element(x,y)==TRUE)]
              return(x1)
              },spw_matrix_rnames)

   retain_drug<-NULL
    for(i in 1:length(Target_upregulation)){
      if((length(Target_upregulation[[i]])>=drug.min.sz)&(length(Target_upregulation[[i]])<=drug.max.sz)){
        if((length(Target_downregulation[[i]])>=drug.min.sz)&(length(Target_downregulation[[i]])<=drug.max.sz)){
          retain_drug<-c(retain_drug,i)
        }
      }
    }
  Target_upregulation<-Target_upregulation[retain_drug]
  Target_downregulation<-Target_downregulation[retain_drug]
  result<-list(Target_upregulation=Target_upregulation,Target_downregulation=Target_downregulation)
  return(result)
}
