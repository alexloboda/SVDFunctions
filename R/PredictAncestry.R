#' SVD-based ancestry prediction
#' 
#' Requires several left singular vector bases from every ancestry to be detected. 
#' 1000 genomes data could be used as training set for creation of reference SVD bases. 
#' Sample genotypes vector is reconstructed from every supplied basis and basis
#'  with a smallest relative residual vector norm is chosen as sample’s ancestry.
#' @param gmatrix Vector of genotypes for a sample
#' @param ReferenceUList List object containing matrices of the left singular
#'  vectors for every ancestry reference
#' @param SV Number of singular vectors to use
#' @param AncestryList Vector of ancestry names
#' @export
PredictAncestry <- function(gmatrix, ReferenceUList, SV, AncestryList){
  if(class(ReferenceUList)!="list"){
    stop("Collection of U-bases must be supplied as list object")
  }
  if(length(which(is.na(gmatrix)==TRUE,arr.ind = TRUE))>0){
    stop("Missing values in genotype matrix. Use ReplaceMissing() function to fix")
  }
  if(length(ReferenceUList)<2){
    stop("More than 1 basis must be supplied")
  }
  if(length(ReferenceUList)!=length(AncestryList)){
    stop("Different length of bases list and ancestries vector")
  }
  resid<-vector("list",length(ReferenceUList))
  for(i in 1:length(ReferenceUList)){
    resid[[i]]<-ParallelResidEstimate(gmatrix = gmatrix,svdReference = ReferenceUList[[i]],nSV = SV)
  }
  resid<-do.call(rbind,resid)
  rownames(resid)<-AncestryList
  svd.pred.anc<-c()
  for(i in 1:ncol(resid)){
    svd.pred.anc<-c(svd.pred.anc,AncestryList[which.min(resid[,i])])
  }
  names(svd.pred.anc)<-colnames(gmatrix)
  return(svd.pred.anc)
}
