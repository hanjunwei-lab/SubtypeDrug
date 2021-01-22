##' PrioSubtypeDrug
##'
##' @title Prioritization of candidate cancer subtype-specific drugs (PrioSubtypeDrug)
##' @description Integrating drug, gene, and subpathway data to identify drugs specific to cancer subtypes.
##' @param expr Matrix of gene expression values (rows are genes, columns are samples).
##' @param input.cls Input sample subtype class vector file in CLS format.
##' @param control.label In the CLS file of `input.cls`, the label of the control sample.
##' @param subpathway.list A list.  The subpathway list data is mined from KEGG data is stored in the package `SubtypeDrugData` and can be downloaded through the connection \url{https://github.com/hanjunwei-lab/SubtypeDrugData}.
##'  The gene tags included in the subpathway list data should be consistent with those in the gene expression profile. The package `SubtypeDrugData` provides two choices that include the Entrezid and Symbol tags of the gene.
##'  Users can also enter their own pathway or gene set list data.
##' @param spw.min.sz Removes subpathways that contain fewer genes than `spw.min.sz` (default: 10).
##' @param spw.max.sz Removes subpathways that contain more genes than `spw.max.sz` (default: Inf).
##' @param spw.score.method Method to employ in the estimation of subpathway
##' enrichment scores per sample. By default this is set to `gsva` (Hänzelmann
##' et al, 2013) and other options are `ssgsea` (Barbie et al, 2009).
##' @param kcdf Character string denoting the kernel to use during the
##' non-parametric estimation of the cumulative distribution function of
##' expression levels across samples when `spw.score.method="gsva"`. By default,
##' `kcdf="Gaussian"` which is suitable when input expression values are
##' continuous, such as microarray fluorescent units in logarithmic scale,
##' RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are
##' integer counts, such as those derived from RNA-seq experiments, then this
##' argument should be set to `kcdf="Poisson"`.
##' @param drug.spw.data A list data of drug regulation. The drug subpathway association data we constructed is stored in package `SubtypeDrugData` and can be downloaded via connection \url{https://github.com/hanjunwei-lab/SubtypeDrugData}.
##' If the input is user-defined drug regulation data, the data should be a list data with each drug as its element. Each drug also contains `Target_upregulation` and `Target_downregulation` subpathway or gene set.
##' Subpathway or gene set contained in drug regulation data should exist in input data of parameter `subpathway.list`.
##' @param drug.spw.p.val.th Parameter used only when `drug.spw.data="DrugSpwData"`. According to the threshold of the significant P value
##' set by parameter `drug.spw.p.val.th` (default: 0.05), the drug up-regulation and down-regulatory subpathways were screened.
##' @param drug.spw.min.sz A numeric. The drug regulated subpathways intersects with the subpathways in the
##' subpathway activity profile. Then drugs with less than `drug.spw.min.sz` (default: 10) up- or down-regulated subpathways are removed.
##' @param drug.spw.max.sz A numeric. Similar to parameter `drug.spw.min.sz`, drugs with more
##' than `drug.spw.max.sz` (default: Inf) up- or down-regulated subpathways are removed.
##' @param weighted.drug.score A boolean values determines the method for
##' calculating the normalized drug-disease reverse association score of the drug for each sample.
##' `weighted.drug.score=TRUE` (default): KS random walk statistic with individualized subpathway activity aberrance score as weight was used to
##' calculate the normalized drug-disease reverse association score.
##' `weighted.drug.score=FALSE`: Similar to `CMap` (Lamb et al., 2006), no weight is needed, and the normalized drug-disease
##' reverse association score is calculated by the rank of the individualized subpathway activity aberrance score.
##' @param nperm Number of random permutations (default: 1000).
##' @param parallel.sz Number of processors to use when doing the calculations in
##' parallel (default value: 1). If parallel.sz=0,
##' then it will use all available core processors unless we set this argument
##' with a smaller number.
##' @param E_FDR Significance threshold for E_FDR for drugs (default: 0.05)
##' @param S_FDR Significance threshold for S_FDR for drugs (default: 0.001)
##' @details
##' First, the function PrioSubtypeDrug uses the `GSVA` or `ssgsea` method to
##' convert the disease gene expression profile into subpathway activity profile. Parameters `subpathway.list`, `spw.min.sz` and `spw.max.sz` are used
##' to process the subpathway list data. `spw.score.method` and `kcdf` are used to control the method of constructing the subpathway activity score
##' profile.
##'   Individualized subpathway activity aberrance score was estimated using the mean and standard deviation of the Control samples.
##' Subpathways of each cancer sample are ordered in a ranked list according to individualized subpathway activity aberrance score.
##' Next, we calculate the normalized drug-disease reverse association score by enriching drug regulated subpathway tags to the subpathway ranked list.
##' Finlly, all drug-regulated subpathways are enriched into each cancer sample to obtain a normalized drug-disease reverse association score matrix.
##' The `drug.spw.p.val.th`, `drug.spw.min.sz` and `drug.spw.max.sz` is used to screen the drug regulated subpathway set.
##'   If user-defined drug targeting data is used, drug regulated `Target_upregulation` and `Target_downregulation` should already be defined in the data.
##' The `weighted.drug.score` to control the method of calculating the normalized drug-disease reverse association score.
##'   Finally, empirical sample-based permutation test procedure to obtain significative cancer subtype specific drugs.
##' For samples containing only cancer and Control, the subpathways are ranked according to the difference in activity between cancer and Control samples.
##' Subsequently, the subpathway set of drug up- and down-regulated is enriched to the ranking list of subpathway to evaluate the normalized drug-disease reverse association score and
##' subpathway-based permutation test procedure to calculate significance.
##' The subpathway list data and drug subpathway associated data set is stored in package `SubtypeDrugData` and
##' can be obtained on \url{https://github.com/hanjunwei-lab/SubtypeDrugData}.
##' @return A list contains the result table of drug scoring and significance, a subpathway activity score matrix, a normalized drug-disease reverse
##' association score matrix, sample information, and user set parameter information.
##' @author Xudong Han,
##' Junwei Han,
##' Chonghui Liu
##' @examples
##' require(GSVA)
##' require(parallel)
##' ## Get simulated breast cancer gene expression profile data.
##' Geneexp<-get("Geneexp")
##' ## Obtain sample subtype data and calculate breast cancer subtype-specific drugs.
##' \donttest{Subtype<-system.file("extdata", "Subtype_labels.cls", package = "SubtypeDrug")}
##'
##' ## Subpathway list data and drug subpathway association data
##' ## were stored in packet `SubtypeDrugData`.
##' ## `SubtypeDrugData` has been uploaded to the github repository.
##' ## If subpathway list data and drug subpathway association data are needed,
##' ## users can download and install through `install_github` function and
##' ## set parameter url=""hanjunwei-lab/SubtypeDrugData".
##' ## After installing and loading package `SubtypeDrugData`,
##' ## users can use the following command to get the data.
##' ## Get subpathway list data.
##' ## If the gene expression profile contains gene Symbol.
##' \donttest{data(SpwSymbolList)}
##' ## If the gene expression profile contains gene Entrezid.
##' \donttest{data(SpwEntrezidList)}
##' ## Get drug subpathway association data.
##' \donttest{data(DrugSpwData)}
##'
##' ## Identify breast subtype-specific drugs.
##' \donttest{Subtype_drugs<-PrioSubtypeDrug(Geneexp,Subtype,"Control",SpwSymbolList,drug.spw.data=DrugSpwData,
##'                                          E_FDR=1,S_FDR=1)}
##'
##' ## Identify breast cancer-related drugs in only two types of samples: breast cancer and control.
##' Cancer<-system.file("extdata", "Cancer_normal_labels.cls", package = "SubtypeDrug")
##' \donttest{Disease_drugs<-PrioSubtypeDrug(Geneexp,Cancer,"Control",SpwSymbolList,drug.spw.data=DrugSpwData,
##'                                          E_FDR=1,S_FDR=1)}
##'
##' ## The function PrioSubtypeDrug() can also support user-defined data.
##' Geneexp<-get("GeneexpT")
##' ## User-defined drug regulation data should resemble the structure below
##' UserDS<-get("UserDST")
##' str(UserDS)
##' ## Need to load gene set data consistent with drug regulation data.
##' UserGS<-get("UserGST")
##' str(UserGS)
##' Drugs<-PrioSubtypeDrug(Geneexp,Cancer,"Control",UserGS,spw.min.sz=1,
##'                        drug.spw.data=UserDS,drug.spw.min.sz=1,
##'                        nperm=10,E_FDR=1,S_FDR=1)
##' @importFrom parallel parLapply
##' @importFrom parallel detectCores
##' @importFrom parallel makeCluster
##' @importFrom parallel clusterExport
##' @importFrom parallel stopCluster
##' @importFrom GSVA gsva
##' @importFrom stats p.adjust
##' @importFrom stats pnorm
##' @export

PrioSubtypeDrug<-function(expr,input.cls="",control.label="",subpathway.list,
                spw.min.sz=10,spw.max.sz=Inf,spw.score.method="gsva",kcdf="Gaussian",
                drug.spw.data,drug.spw.p.val.th=0.05,
                drug.spw.min.sz=10,drug.spw.max.sz=Inf,
                weighted.drug.score=TRUE,nperm=1000,parallel.sz=1,E_FDR=0.05,S_FDR=0.001){

  haveGSVA <- isPackageLoaded("GSVA")
  if(haveGSVA==FALSE){
    stop("The 'GSVA' library, should be loaded first")
  }
  haveParallel <- isPackageLoaded("parallel")
  if(haveParallel==FALSE){
    stop("The 'parallel' library, should be loaded first")
  }

  spw_matrix_y<-gsva(expr,subpathway.list,method=spw.score.method,kcdf=kcdf,
                     min.sz=spw.min.sz,max.sz=spw.max.sz,
                     parallel.sz=parallel.sz,verbose=FALSE)
  spw_matrix_rnames<-row.names(spw_matrix_y)
  sampleid<-colnames(spw_matrix_y)

  if(is.list(input.cls)) {
    CLS <- input.cls
  } else {
    CLS <- ReadClsFile(file=input.cls)
  }
  phen <- CLS$phen
  samples.v1 <-CLS$class.labes
  SmaplePhenotype<-data.frame(sampleId=sampleid,sampleSubtype=samples.v1,stringsAsFactors =FALSE)
  control_index<-which(samples.v1==control.label)
  phen<-phen[phen!=control.label]
  samples.v<-samples.v1[-control_index]

  drug_target_data<-getDrugSpw(drug.spw.data,spw_matrix_rnames,drug.spw.p.val.th,drug.spw.min.sz,drug.spw.max.sz)

  up_signature1<-sapply(drug_target_data$Target_upregulation, function(x){
        pj<-paste(x,collapse =" ")
      })
  down_signature1<-sapply(drug_target_data$Target_downregulation, function(x){
        pj<-paste(x,collapse =" ")
      })
  Parameters<-list(control.label=control.label,spw.min.sz=spw.min.sz,spw.max.sz=spw.max.sz,spw.score.method=spw.score.method,kcdf=kcdf,
                    drug.spw.p.val.th=drug.spw.p.val.th,drug.spw.min.sz=drug.spw.min.sz,drug.spw.max.sz=drug.spw.max.sz,weighted.drug.score=weighted.drug.score,nperm=nperm,parallel.sz=parallel.sz)

   if(length(phen)==1){
      fc<-apply(spw_matrix_y, 1, function(x){
        foldc<-mean(x[which(samples.v1!=control.label)])-mean(x[control_index])
        return(foldc)
      })
      drugs<-getDrugMatrix(fc,drug_target_data,weighted.drug.score)
      drugs1<-drugs
      zdz<-max(drugs)
      zxz<-min(drugs)
      for(i in 1:length(drugs)){
       if(drugs[i]>0){
         drugs[i]<-drugs[i]/abs(zdz)
       }else{
          drugs[i]<-drugs[i]/abs(zxz)
        }
     }

      if(parallel.sz==1){
        fc1<-fc
        rdmean_matrix<-NULL
        for(n in 1:nperm){
          rdsamples<-sample(c(1:length(fc1)),size = length(fc1),replace = FALSE)
          names(fc1)<-names(fc1)[rdsamples]
          rdmatrix<-getDrugMatrix(fc1,drug_target_data,weighted.drug.score)
          rdmean_matrix<-cbind(rdmean_matrix,rdmatrix)
        }
      }else{
        if(parallel.sz==0){
         nCores <- parallel::detectCores()
        }else{
         nCores <- parallel.sz
      }
      cl <- makeCluster(nCores)
      clusterExport(cl,c("CalculateSES","getDrugMatrix"))
      rdmean_matrix<-parLapply(cl,1:nperm,function(n,fc1,weighted.drug.score,drug_target_data){
        rdsamples<-sample(c(1:length(fc1)),size = length(fc1),replace = FALSE)
        names(fc1)<-names(fc1)[rdsamples]
        rdmatrix<-getDrugMatrix(fc1,drug_target_data,weighted.drug.score)
        return(rdmatrix)
      },fc,weighted.drug.score,drug_target_data)
      stopCluster(cl)
      rdmean_matrix<-do.call("cbind",rdmean_matrix)
      }

      p_values<-sapply(c(1:length(drugs1)), function(i){
        s_t<-abs(drugs1[i])
        p.val<-pnorm(-s_t,mean = mean(rdmean_matrix[i,]), sd = sd(rdmean_matrix[i,]),lower.tail = TRUE)+pnorm(s_t,mean = mean(rdmean_matrix[i,]), sd = sd(rdmean_matrix[i,]),lower.tail = F)
        return(p.val)
      })

      fdr<-p.adjust(p_values,"BH",length(p_values))

      result<-data.frame(Drug=names(drugs),
                          Target_upregulation=up_signature1,
                          Target_downregulation=down_signature1,
                          SDS=signif(drugs,digits=3),E_Pvalue=signif(p_values,digits=3),E_FDR=signif(fdr,digits=3),stringsAsFactors=FALSE)
      result<-result[which(result[,6]<=E_FDR),]
      result<-list(result,spw_matrix_y,SmaplePhenotype,Parameters)
      names(result)<-c(phen,"SubpathwayMatrix","SampleInformation","Parameter")
      return(result)
   }else{

     sy_spw_matrix<-AccumulateNormal(spw_matrix_y,control_index)

     drug_sample_matrix<-getDrugMatrix(sy_spw_matrix,drug_target_data,weighted.drug.score)
     drug_sample_matrix<-t(apply(drug_sample_matrix,1,function(x){
      zdz<-max(x)
      zxz<-min(x)
      for(i in 1:length(x)){
       if(x[i]>0){
         x[i]<-x[i]/abs(zdz)
       }else{
          x[i]<-x[i]/abs(zxz)
        }
      }
       return(x)
     }))

     pn<-length(phen)
     drug_true_s<-sapply(1:pn, function(p){
       phenindex<-which(samples.v==phen[p])
       DCS<-apply(drug_sample_matrix[,phenindex], 1,mean)
       return(DCS)
     })


     if(parallel.sz==0){
       nCores <- parallel::detectCores()
     }else{
       nCores <- parallel.sz
     }
     cl <- makeCluster(nCores)
     rdmean_matrix<-parLapply(cl,1:nperm,function(n,drug_matrix,samples.v,phen){
        rdsamples<-sample(samples.v,size = length(samples.v),replace = FALSE)
        rdmeans<-NULL
        for(i in 1:length(phen)){
          ppindex<-which(rdsamples==phen[i])
          rdmeans<-cbind(rdmeans,apply(drug_matrix[,ppindex], 1, mean))
        }
        return(rdmeans)
      },drug_sample_matrix,samples.v,phen)

      rdmean_matrix<-do.call("cbind",rdmean_matrix)

      sub_fc<-NULL
      sub_score<-NULL
      for(p in 1:length(phen)){
        fc<-apply(spw_matrix_y, 1, function(x){
          foldc<-mean(x[which(samples.v1==phen[p])])-mean(x[which(samples.v1==control.label)])
          return(foldc)
       })
        sub_score<-cbind(sub_score,getDrugMatrix(fc,drug_target_data,weighted.drug.score))
        sub_fc<-cbind(sub_fc,fc)
      }

      clusterExport(cl,c("CalculateSES","getDrugMatrix"))
      subrd_matrix<-parLapply(cl,1:nperm,function(n,fc1,weighted.drug.score,drug_target_data){
        row.names(fc1)<-sample(row.names(fc1),size = nrow(fc1),replace = FALSE)
        rdmatrix<-apply(fc1,2,function(x){
          return(getDrugMatrix(x,drug_target_data,weighted.drug.score))
        })
        return(rdmatrix)
      },sub_fc,weighted.drug.score,drug_target_data)
      subrd_matrix<-do.call("cbind",subrd_matrix)
      stopCluster(cl)

      sub_p<-NULL
      sub_fdr<-NULL
      for(p in 1:length(phen)){
        if(p==length(phen)){
          sx<-which(c(1:ncol(subrd_matrix))%%length(phen)==0)
        }else{
          sx<-which(c(1:ncol(subrd_matrix))%%length(phen)==p)
        }
        subrd_matrix1<-subrd_matrix[,sx]
        sub_p1<-sapply(c(1:nrow(sub_score)), function(i){
          s_t<-abs(sub_score[i,p])
          p.val<-pnorm(-s_t,mean = mean(subrd_matrix1[i,]), sd = sd(subrd_matrix1[i,]),lower.tail = TRUE)+pnorm(s_t,mean = mean(subrd_matrix1[i,]), sd = sd(subrd_matrix1[i,]),lower.tail = F)
          return(p.val)
        })
         sub_fdr1<-p.adjust(sub_p1,"BH",length(sub_p1))
         sub_p<-cbind(sub_p,sub_p1)
         sub_fdr<-cbind(sub_fdr,sub_fdr1)
      }

      colnames(drug_true_s)<-phen
      p_values<-sapply(c(1:nrow(drug_true_s)),function(d){
        p_val_v<-NULL
        for(i in 1:pn){
          s_t<-abs(drug_true_s[d,i])
          p.val<-pnorm(-s_t,mean = mean(rdmean_matrix[d,]), sd = sd(rdmean_matrix[d,]),lower.tail = TRUE)+pnorm(s_t,mean = mean(rdmean_matrix[d,]), sd = sd(rdmean_matrix[d,]),lower.tail = F)
          p_val_v<-c(p_val_v,p.val)
        }
        return(p_val_v)
     })

      result<-list()
      for(i in 1:pn){
        drug_true_s1<-signif(drug_true_s[,i],digits=3)
        fdr<-p.adjust(p_values[i,],"BH",ncol(p_values))
        result1<-data.frame(Drug=row.names(drug_sample_matrix),
                         Target_upregulation=up_signature1,
                         Target_downregulation=down_signature1,
                         SDS=drug_true_s1,E_Pvalue=sub_p[,i],E_FDR=sub_fdr[,i],S_Pvalue=signif(p_values[i,],digits=3),S_FDR=signif(fdr,digits=3),stringsAsFactors=FALSE)
        result[[i]]<-result1[which(result1[,6]<=E_FDR&result1[,8]<=S_FDR),]
      }

    result<-c(result,list(drug_sample_matrix,spw_matrix_y,SmaplePhenotype,Parameters))
    names(result)<-c(phen,"DrugMatrix","SubpathwayMatrix","SampleInformation","Parameter")
    return(result)
   }
}
