#' Plot details of a performed matching.
#' 
#' The function is aimed at plotting the both cohorts (case and control) 
#' in terms of PCA, the MVN distribution is shown as a set of loadings 
#' projected into subspace spanned by PC loadings of control cohort.
#' @param popObj a subpopulation from instance read by function readInstanceFromYml.
#' For an instance object inst an objects like inst$population[[1]] can be passed.
#' @param variants list of variants gathered from instance object.
#' @param controlsU a numerical matrix carrying vectors of U component from SVD.
#' It may contain superset of variants provided by variants parameter. 
#' Rownames should be set to variant signatures.
#' @param meanControl mean genotype value of variants in the order of rows of controlsU matrix.
#' @param casesGenotypeMatrix case genotype matrix, row names must be set to variant signatures.
#' @param controlsGenotypeMatrix control genotype matrix, row names must be set to variant signatures.
#' colnames must be set to sample names.
#' @param matchedControls names of matched controls.
#'
#' @export
plotMatching <- function(popObj, variants, controlsU, meanControl, 
                         casesGenotypeMatrix, controlsGenotypeMatrix, 
                         matchedControls) {
  reduced_vars <- intersect(variants, rownames(controlsU))
  names(meanControl) <- rownames(controlsU)
  controlsU <- controlsU[reduced_vars, ]
  meanControl <- meanControl[reduced_vars]
  
  casesGenotypeMatrix <- casesGenotypeMatrix[reduced_vars, ]
  controlsGenotypeMatrix <- controlsGenotypeMatrix[reduced_vars, ]
  
  Uinv = pracma::pinv(controlsU)
  casesGenotypeMatrix <- casesGenotypeMatrix - meanControl
  controlsGenotypeMatrix <- controlsGenotypeMatrix - meanControl
  
  rsCases <- Uinv %*% casesGenotypeMatrix
  rsControls <- Uinv %*% controlsGenotypeMatrix
  
  casesMean <- popObj$mean 
  casesLoadings <- popObj$US
 
  function(pc1, pc2) {
    casePCLs <- casesLoadings[c(pc1, pc2), ]
    arrows <- data.frame(startx = casesMean[pc1], 
                         starty = casesMean[pc2], 
                         endx = casesMean[pc1] + casePCLs[1, ], 
                         endy = casesMean[pc2] + casePCLs[2, ])
    
    covPCL <- t(casesLoadings) %*% casesLoadings
    effect <- expm::sqrtm(covPCL)
    main_arrows <- data.frame(startx = casesMean[pc1], 
                              starty = casesMean[pc2], 
                              endx = casesMean[pc1] + effect[1,], 
                              endy = casesMean[pc2] + effect[2,])
    
    df <- data.frame(PC1 = rsControls[pc1,], PC2 = rsControls[pc2,], color = "control")
    rownames(df) <- colnames(controlsGenotypeMatrix)
    df[matchedControls, ]$color <- "selected"
    
    casesDF <- data.frame(PC1 = rsCases[pc1, ], PC2 = rsCases[pc2, ], color = "case")
    
    ggplot() + 
     geom_point(aes(PC1, PC2, colour = "matched controls"), df[df$color == "selected",], alpha = 0.1, color = "darkred") +
     geom_point(aes(PC1, PC2, colour = "available controls"), df[df$color == "control",], alpha = 0.01)  +
     geom_point(aes(PC1, PC2, colour = "cases"), casesDF, alpha = 1.1, color = "darkblue") + 
     geom_segment(aes(x = startx, y = starty, xend = endx, yend = endy, color = "PC loadings projections"), arrows, 
                  arrow = arrow(length = unit(0.1, "inches")), size = 1.5, alpha = 0.8) +
     geom_segment(aes(x = startx, y = starty, xend = endx, yend = endy, color = "novel PC loadings"), main_arrows, 
                  arrow = arrow(length = unit(0.1, "inches")), size = 1.5) + 
     xlab(paste0("PC", pc1)) + ylab(paste0("PC", pc2)) + theme_bw() 
  } 
}