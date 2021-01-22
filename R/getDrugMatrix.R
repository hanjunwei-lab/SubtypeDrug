##' getDrugMatrix
##'
##'
##' @title SubtypeDrug internal function
##' @description Obtaining drug-disease reverse association score matrix.
##' @param spw_matrix A subpathway activity profile. rows are subpathwyas, columns are samples.
##' @param drug_target_data A list. A list stores a collection of drug up- and down-regulated subpathways.
##' @param weighted.score A binary value of 0 or 1. If the `weighted.score` = 1, the drug reverse association score will be weighted by the subpathway activity.
##' @return A matrix.
##' @examples
##' \donttest{require(GSVA)}
##' \donttest{Geneexp<-get("GeneexpT")}
##' \donttest{UserDS<-get("UserDST")}
##' \donttest{UserGS<-get("UserGST")}
##' \donttest{spw_matrix<-gsva(Geneexp,UserGS,verbose=F)}
##' \donttest{x<-getDrugMatrix(spw_matrix,UserDS,weighted.score=1)}
##' @author Xudong Han,
##' Junwei Han,
##' Chonghui Liu
##' @export
getDrugMatrix<-function(spw_matrix,drug_target_data,weighted.score){
  ism<-is.matrix(spw_matrix)
  if(ism==TRUE){
    if(weighted.score==FALSE){
       drug.sample.matrix <- apply(spw_matrix, 2, function(x){
          fc_final_sort <- sort(x, decreasing = TRUE)
          ks_up_list <- lapply(drug_target_data$Target_upregulation, function(y) {
             n <- length(fc_final_sort)
            t <- length(y)
            V <- match(y, names(fc_final_sort))
            V <- sort(V)

             j <- 1:t
             a <- j/t - V/n
             a <- max(a)


             b <- V/n - (j - 1)/t
             b <- max(b)

             if (a >= b) {
              ks_up = a
             } else {
               ks_up = -b
            }
             return(ks_up)
           })

          ks_down_list <- lapply(drug_target_data$Target_downregulation, function(y) {
             n <- length(fc_final_sort)
             t <- length(y)
             V <- match(y, names(fc_final_sort))
             V <- sort(V)

             j <- 1:t
             a <- j/t - V/n
             a <- max(a)

             b <- V/n - (j - 1)/t
             b <- max(b)

             if (a >= b) {
               ks_down = a
            } else {
               ks_down = -b
            }
            return(ks_down)
          })

           s <- unlist(ks_up_list) - unlist(ks_down_list)
           names(s) <- names(ks_up_list)
           return(s)
        })
    }else{
      drug.sample.matrix<- apply(spw_matrix, 2, function(x) {
        x <- sort(x, decreasing = T)
        spw_matrix_rnames_rank<-names(x)
        up.drug.sample.score <- sapply(drug_target_data$Target_upregulation, function(z) {
          tag.indicator <- sign(match(spw_matrix_rnames_rank, z, nomatch = 0))
          up.ES <- CalculateSES(tag.indicator,correl.vector = x)
          return(up.ES)
        })
        down.drug.sample.score <- sapply(drug_target_data$Target_downregulation, function(z) {
          tag.indicator <- sign(match(spw_matrix_rnames_rank, z, nomatch = 0))
          down.ES <- CalculateSES(tag.indicator,correl.vector = x)
          return(down.ES)
        })
        drug.sample.score <- up.drug.sample.score - down.drug.sample.score
        return(drug.sample.score)
      })
    }
    return(drug.sample.matrix)
  }else{
    if(weighted.score==FALSE){
      fc_final_sort <- sort(spw_matrix, decreasing = TRUE)
      ks_up_list<-lapply(drug_target_data$Target_upregulation, function(y){
         n <- length(fc_final_sort)
         t <- length(y)
         V <- match(y, names(fc_final_sort))
         V <- sort(V)

         j <- 1:t
         a <- j/t - V/n
         a <- max(a)


         b <- V/n - (j - 1)/t
         b <- max(b)

         if (a >= b) {
           ks_up = a
         } else {
           ks_up = -b
         }
         return(ks_up)
       })

      ks_down_list <- lapply(drug_target_data$Target_downregulation, function(y) {
         n <- length(fc_final_sort)
         t <- length(y)
         V <- match(y, names(fc_final_sort))
         V <- sort(V)

         j <- 1:t
         a <- j/t - V/n
         a <- max(a)

         b <- V/n - (j - 1)/t
         b <- max(b)

         if (a >= b) {
           ks_down = a
         } else {
           ks_down = -b
         }
         return(ks_down)
       })
      drug.sample.matrix <- unlist(ks_up_list) - unlist(ks_down_list)
      names(drug.sample.matrix) <- names(ks_up_list)

    }else{
      x <- sort(spw_matrix, decreasing = T)
      spw_matrix_rnames_rank<-names(x)
      up.drug.sample.score <- sapply(drug_target_data$Target_upregulation, function(z) {
        tag.indicator <- sign(match(spw_matrix_rnames_rank, z, nomatch = 0))
        up.ES <- CalculateSES(tag.indicator,correl.vector = x)
        return(up.ES)
      })
      down.drug.sample.score <- sapply(drug_target_data$Target_downregulation, function(z) {
        tag.indicator <- sign(match(spw_matrix_rnames_rank, z, nomatch = 0))
        down.ES <- CalculateSES(tag.indicator,correl.vector = x)
        return(down.ES)
      })
      drug.sample.matrix <- up.drug.sample.score - down.drug.sample.score
    }
    return(drug.sample.matrix)
  }


}
