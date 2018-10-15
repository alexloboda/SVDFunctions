#' Replace of missing values in genotype matrix
#' 
#' SVD does not tolerated missing values in a matrix. This function replaces 
#' missing entries in genotype matrix with average genotype for each DNA 
#' variant rounded to an integer.
#' @param GenotypeMatrix Matrix of genotypes (class – Matrix)
#' @export
ReplaceMissing <- function (GenotypeMatrix) 
{
    k <- which(is.na(GenotypeMatrix), arr.ind = T)
    if (length(k) > 0) {
        GenotypeMatrix[k] <- round(rowMeans(GenotypeMatrix, na.rm = T)[k[, 1]])
    }
    return(GenotypeMatrix)
}
