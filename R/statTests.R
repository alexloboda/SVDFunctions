#' Weighted Linear Regression for Genotype-Phenotype Association
#'
#' Performs a fast linear regression of phenotype on genotype using provided case/control counts.
#'
#' @param case_0 Number of cases with genotype 0
#' @param case_1 Number of cases with genotype 1
#' @param case_2 Number of cases with genotype 2
#' @param control_0 Number of controls with genotype 0
#' @param control_1 Number of controls with genotype 1
#' @param control_2 Number of controls with genotype 2
#'
#' @return A list with elements \code{beta} (regression coefficient) and \code{pval} (p-value)
#' @examples
#' genoPhenoLM(10, 20, 30, 15, 25, 35)
#' @export
genoPhenoLM <- function(case_0, case_1, case_2, control_0, control_1, control_2) {
  # Argument checks
  counts <- c(case_0, case_1, case_2, control_0, control_1, control_2)
  if (any(is.na(counts)) || any(counts < 0)) {
    stop("All counts must be non-negative and not NA.")
  }
  if (sum(counts) == 0) {
    stop("At least one count must be positive.")
  }
  # Prepare data
  geno <- rep(0:2, 2)
  phe <- c(rep(0, 3), rep(1, 3))
  ws <- c(control_0, control_1, control_2, case_0, case_1, case_2)
  table <- data.frame(geno = geno, phe = phe)
  # Fit model
  m <- lm(phe ~ geno, data = table, weights = ws)
  m$df.residual <- sum(ws) - 2
  ks <- summary(m)$coefficients
  list(
    beta = ks["geno", "Estimate"],
    pval = ks["geno", "Pr(>|t|)"]
  )
}

#' Calculate Genomic Inflation Factor
#'
#' Calculates the genomic inflation factor (lambda) from a vector of p-values.
#' If the input p-values are zeros or ones, they are
#'  replaced with uniform random values
#' between zero and the first non-zero p-value, and between
#' the last non-one p-value and one, respectively.
#'
#' The fraction of smallest p-values to consider
#'  for lambda calculation can be adjusted
#' This may be useful in the case of applying to results of
#' Fisher's exact test, where the p-values are not
#' ideally uniformly distributed.
#' @param pvals Numeric vector of p-values
#' @param fractionSmallest fracrtion of smallest p-values 
#' to consider for lambda calculation
#' @return Numeric value of the genomic inflation factor
#' @examples
#' pvals <- c(0.01, 0.05, 0.1, 0.2)
#' lambda_val <- calcLambda(pvals)
#' @export
calcLambda <- function(pvals, fractionSmallest = 1.0) {
  if (length(pvals) == 0 || all(is.na(pvals))) {
    return(NA)
  }
  pvals <- pvals[!is.na(pvals)]
  # pvals zero and one replace by unform between zero and the first non-zero, the last non-one and one correspondingly
  EPS <- 1e-20
  mask_zero <- pvals < EPS
  mask_one <- pvals > (1 - EPS)
  pvals[mask_zero] <- runif(sum(mask_zero), 0, min(pvals[pvals >= EPS], na.rm = TRUE))
  pvals[mask_one] <- runif(sum(mask_one), max(pvals[pvals <= (1 - EPS)], na.rm = TRUE), 1)

  n_obs <- min(length(pvals), max(1, round(length(pvals) * fractionSmallest)))
  pvals <- sort(pvals)[1:n_obs]

  if (length(pvals) == 0) {
    return(NA)
  }
  chisq <- qchisq(1 - pvals, df = 1)
  expected_median <- qchisq(1.0 - fractionSmallest / 2.0, df = 1)
  lambda_val <- median(chisq) / expected_median
  if (is.nan(lambda_val) || is.infinite(lambda_val)) {
    return(NA)
  }
  return(lambda_val)
}

#' QQ Plot for P-values
#' 
#' Generates a QQ plot for the provided p-values.
#' @param pvals Numeric vector of p-values
#' @param title Title for the plot
#' @return A ggplot object
#' @export
qqPlot <- function(pvals, title = "QQ Plot") {
  lambda_val <- calcLambda(pvals)

  pvals1 <- pvals
  pvals2 <- ppoints(length(pvals))
  pvals.data <- data.frame(X = -log10(sort(pvals)),
                           Y = -log10(sort(pvals2)), 
                           expected = -log10(ppoints(length(pvals1))),
                           clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:length(pvals1),
                                                   shape2 = length(pvals1):1)),
                           cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:length(pvals1),
                                                   shape2 = length(pvals1):1)))

  ggplot2::ggplot(pvals.data) +
  ggplot2::ggtitle(title) + 
  ggplot2::geom_abline(intercept = 0,
              slope = 1, 
              alpha = 0.5,color=adjustcolor("grey",alpha.f = 1)) +
  ggplot2::geom_line(ggplot2::aes(expected, cupper), linetype = 2,color=adjustcolor("grey",alpha.f = 0.5)) +
  ggplot2::geom_line(ggplot2::aes(expected, clower), linetype = 2,color=adjustcolor("grey",alpha.f = 0.5)) +
  ggplot2::geom_point(data = pvals.data, ggplot2::aes(x = expected, 
                                    y = X,
                                    colour="External"),
             size=2,shape=20,stroke=0) +
  ggplot2::xlab(expression(paste("Expected -log"[10], plain(P)))) +
  ggplot2::ylab(expression(paste("Observed -log"[10], plain(P))))+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.text=ggplot2::element_text(size=8),
        plot.title=ggplot2::element_text(family="Palatino", size=8),
        axis.title=ggplot2::element_text(size=9),
        legend.text=ggplot2::element_text(size=8),
        legend.title=ggplot2::element_text(size=7),
        plot.margin = ggplot2::unit(c(.1,.3,.1,.1), "cm"),
        aspect.ratio = 1,
        legend.key.size = ggplot2::unit(0.5,"cm"),
        legend.position = c(0.35,0.8),
        legend.background = ggplot2::element_rect(color = "black", fill = "white", size = 0.2, linetype = "solid"),
        legend.margin = ggplot2::margin(t = .1,r = .1,b = .1,l = .1,unit='cm'))+
  ggplot2::labs(fill="Genomic inflation")+
  ggplot2::scale_colour_manual(name="Genomic Inflation",
                      guide='legend',
                      labels=c(paste("λ", "=", format(lambda_val, digits = 3))),
                      values=c("black"))
}

#' Vectorized Linear Regression P-value
#'
#' Vectorized version of \code{genoPhenoLM} returning only the p-value.
#'
#' @inheritParams genoPhenoLM
#' @return Numeric p-value
#' @examples
#' # Vectorized usage:
#' genoPhenoVecPval(
#'   case_0 = c(10, 11), case_1 = c(20, 21), case_2 = c(30, 31),
#'   control_0 = c(15, 16), control_1 = c(25, 26), control_2 = c(35, 36)
#' )
#' @export
genoPhenoVecPval <- Vectorize(function(case_0, case_1, case_2, control_0, control_1, control_2) {
  genoPhenoLM(case_0, case_1, case_2, control_0, control_1, control_2)$pval
})
