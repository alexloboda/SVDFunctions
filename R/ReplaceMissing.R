#' Replace of missing values in genotype matrix
#' 
#' SVD does not tolerated missing values in a matrix. This function replaces 
#' missing entries in genotype matrix with average genotype for each DNA 
#' variant rounded to an integer.
#' @param genotypeMatrix Matrix of genotypes (class – Matrix)
#' @export
replaceMissing <- function (genotypeMatrix) 
{
  for (i in 1:nrow(genotypeMatrix)) {
    m <- round(mean(genotypeMatrix[i, ], na.rm = TRUE))
    genotypeMatrix[i, is.na(genotypeMatrix[i, ])] <- m
  }
  genotypeMatrix
}
