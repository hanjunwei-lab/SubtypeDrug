##' getDrugStructure
##'
##'
##' @title Get drug chemical structure diagram data
##' @description `getDrugStructure()` outputs the chemical structure graph data of the
##' drug or compound based on the input drug label by the user. The results can be visualized by the `plot` function.
##' @param drug.label A character string of drug label to determine which drug to use for visualization.
##' @param main An overall title for the chemical structure graph.
##' @param sub A sub title for the chemical structure graph.
##' @return A sdfset object.
##' @author Xudong Han,
##' Junwei Han,
##' Chonghui Liu
##' @examples
##' \donttest{require(rvest)}
##' \donttest{require(ChemmineR)}
##' # Plot the chemical structure of drug pirenperone.
##' \donttest{Chem_str<-getDrugStructure(drug.label="pirenperone.")}
##' \donttest{plot(Chem_str)}
##' @importFrom ChemmineR read.SDFset
##' @importFrom rvest html_text
##' @importFrom xml2 read_html
##' @importFrom graphics plot
##' @export

getDrugStructure<-function(drug.label="",main="",sub=""){
  haveChemmineR <- isPackageLoaded("ChemmineR")
  havervest <- isPackageLoaded("rvest")
  if(haveChemmineR==FALSE){
    stop("The 'ChemmineR' library, should be loaded first")
  }
  if(havervest==FALSE){
    stop("The 'rvest' library, should be loaded first")
  }
  drug.label<-unlist(strsplit(drug.label,"\\("))[1]
  drug.label1<-tolower(drug.label)
  Drugs_CID<-get("Drugs_CID")
  drugCid<-Drugs_CID[which(Drugs_CID[,1]==drug.label1),2]
  drug_url<-paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/",drugCid,"/record/SDF/?record_type=2d&response_type=display",sep = "")
  cw<-try(read_html(drug_url))
  if ('try-error' %in% class(cw)){
    stop("Please ensure smooth network connection")
  }
  #drugnr<-read_html(drug_url)
  drugnr<-html_text(cw)
  drugnr<-strsplit(drugnr,"\n")
  drugnr<-unlist(drugnr)
  sdfset <- read.SDFset(drugnr)
  if(main==""){
    sdfset@ID <- drug.label
  }else{
    sdfset@ID <- main
  }
  return(sdfset)
}
