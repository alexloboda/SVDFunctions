#' SVD-based ancestry prediction
#' 
#' Requires several left singular vector bases from every ancestry to be detected. 
#' 1000 genomes data could be used as training set for creation of reference SVD bases. 
#' Sample genotypes vector is reconstructed from every supplied basis and basis
#'  with a smallest relative residual vector norm is chosen as sample’s ancestry.
#' @param genotypeMatrix Vector of genotypes for a sample
#' @param referenceUList List object containing matrices of the left singular
#'  vectors for every ancestry reference
#' @param SV Number of singular vectors to use
#' @param ancestryList Vector of ancestry names
#' @export
predictAncestry <- function(genotypeMatrix, referenceUList, SV, ancestryList){
  gmatrix <- genotypeMatrix
  if(class(referenceUList)!="list"){
    stop("Collection of U-bases must be supplied as list object")
  }
  if(length(which(is.na(gmatrix)==TRUE,arr.ind = TRUE))>0){
    stop("Missing values in genotype matrix. Use ReplaceMissing() function to fix")
  }
  if(length(referenceUList)<2){
    stop("More than 1 basis must be supplied")
  }
  if(length(referenceUList)!=length(ancestryList)){
    stop("Different length of bases list and ancestries vector")
  }
  resid<-vector("list",length(referenceUList))
  for(i in 1:length(referenceUList)){
    resid[[i]]<-ParallelResidEstimate(genotypeMatrix = gmatrix, 
                                      SVDReference = referenceUList[[i]], 
                                      nSV = SV)
  }
  resid<-do.call(rbind,resid)
  rownames(resid)<-ancestryList
  svd.pred.anc<-c()
  for(i in 1:ncol(resid)){
    svd.pred.anc<-c(svd.pred.anc,ancestryList[which.min(resid[,i])])
  }
  names(svd.pred.anc)<-colnames(gmatrix)
  return(svd.pred.anc)
}
