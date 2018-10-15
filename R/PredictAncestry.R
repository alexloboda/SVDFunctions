#' SVD-based ancestry prediction
#' 
#' Requires several left singular vector bases from every ancestry to be detected. 
#' 1000 genomes data could be used as training set for creation of reference SVD bases. 
#' Sample genotypes vector is reconstructed from every supplied basis and basis
#'  with a smallest relative residual vector norm is chosen as sample’s ancestry.
#' @param SampleVector Vector of genotypes for a sample
#' @param ReferenceUList List object containing matrices of the left singular
#'  vectors for every ancestry reference
#' @param SV Number of singular vectors to use
#' @param AncestryList Vector of ancestry names
#' @export
PredictAncestry <- function(SampleVector, ReferenceUList, SV, AncestryList){
  residuals <- rep(0, length(AncestryList))
  for (i in 1:length(AncestryList)){
    residuals[i] <- ComputeResidual.Cross(SampleVector, ReferenceUList[[i]], SV)
  }
  return(AncestryList[which.min(residuals)])
}
