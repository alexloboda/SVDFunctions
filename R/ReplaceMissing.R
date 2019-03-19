#' Replace of missing values in genotype matrix
#' 
#' SVD does not tolerated missing values in a matrix. This function replaces 
#' missing entries in genotype matrix with average genotype for each DNA 
#' variant rounded to an integer.
#' @param genotypeMatrix Matrix of genotypes (class – Matrix)
#' @export
replaceMissing <- function (genotypeMatrix) 
{
    k <- which(is.na(genotypeMatrix), arr.ind = T)
    if (length(k) > 0) {
        genotypeMatrix[k] <- round(rowMeans(genotypeMatrix, na.rm = T)[k[, 1]])
    }
    return(genotypeMatrix)
}
