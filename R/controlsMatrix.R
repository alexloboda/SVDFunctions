#' Write a control matrix to disk.
#'
#' Stores a control genotype matrix in the on-disk format the file backed
#' control selection reads. Rows are variants and every row is stored
#' contiguously, which is what lets \code{\link{prepareControlsSpace}} and
#' \code{\link{matchControlCandidates}} stream the matrix one variant at a time
#' instead of holding it in memory.
#'
#' @param m numeric or integer matrix with variant row names and sample column
#' names.
#' @param path the file to create.
#' @param dtype the on-disk element type. \code{"f64"} stores imputed genotypes
#' exactly. \code{"f32"} halves the file at the cost of single precision, which
#' is far below the resolution of a genotype dosage. \code{"u8"} stores raw
#' genotypes 0, 1, 2 and \code{NA} in one byte each and is the type
#' \code{\link{matchControlCandidates}} reads; it rejects any other value.
#' @return \code{path}, invisibly.
#' @seealso \code{\link{readControlsMatrixInfo}},
#' \code{\link{collapseControlsMatrix}}, \code{\link{selectControlsFromFiles}}
#' @export
writeControlsMatrix <- function(m, path, dtype = c("f64", "f32", "u8")) {
  dtype <- match.arg(dtype)
  stopifnot(is.matrix(m))
  if (!is.numeric(m)) {
    stop("the control matrix must be numeric or integer")
  }
  variants <- rownames(m)
  samples <- colnames(m)
  if (is.null(variants)) {
    stop("the control matrix must have variant row names")
  }
  if (is.null(samples)) {
    stop("the control matrix must have sample column names")
  }
  write_controls_matrix_cpp(m, path.expand(path), dtype, variants, samples)
  invisible(path)
}

#' Read the header of a control matrix file.
#'
#' Reads everything but the data, so it stays cheap no matter how large the
#' matrix is.
#' @param path the control matrix file.
#' @return a list with the element type (\code{dtype}), the dimensions
#' (\code{nrow} variants, \code{ncol} samples or clusters), the
#' \code{rownames} and \code{colnames}, and \code{colweights} -- the number of
#' control samples each column stands for, which is one per column for a
#' genotype matrix and the cluster size for collapsed cluster counts.
#' @seealso \code{\link{writeControlsMatrix}}
#' @export
readControlsMatrixInfo <- function(path) {
  info <- controls_matrix_info_cpp(path.expand(path))
  info$nrow <- as.integer(info$nrow)
  info$ncol <- as.integer(info$ncol)
  info
}

#' Collapse a raw control matrix into per-cluster allele counts.
#'
#' Turns a raw genotype matrix written with \code{dtype = "u8"} into the
#' per-cluster allele counts the matching stage consumes, streaming one variant
#' at a time so neither the input nor the output is ever held in R. The result
#' is both much smaller than the genotype matrix and free of individual
#' genotypes, and it is what \code{\link{matchControlCandidates}} reads fastest.
#'
#' Cluster labels are ordered exactly as \code{\link{prepareControlsSpace}}
#' orders them for the same clustering, so the two stages agree by construction.
#' @param path the source control matrix file, written with \code{dtype = "u8"}.
#' @param outPath the cluster counts file to create.
#' @param controlsClustering cluster names for controls, one per column of the
#' source matrix. When \code{NULL} every sample forms its own cluster. Each
#' cluster must contain at most 255 samples, because a count is stored in one
#' byte.
#' @return invisibly, a list with the output \code{path}, the number of
#' \code{variants}, the number of \code{clusters}, the \code{clusterLabels} and
#' the per-cluster sample sizes (\code{clusterSizes}).
#' @seealso \code{\link{writeControlsMatrix}}, \code{\link{selectControlsFromFiles}}
#' @export
collapseControlsMatrix <- function(path, outPath, controlsClustering = NULL) {
  path <- path.expand(path)
  info <- readControlsMatrixInfo(path)
  if (info$dtype != "u8") {
    stop("collapsing needs a raw genotype matrix written with dtype \"u8\"")
  }
  clusteringInfo <- prepareControlsClustering(controlsClustering, info$colnames,
                                              info$ncol)
  res <- collapse_controls_matrix_cpp(path, clusteringInfo$sampleClusterIds,
                                      clusteringInfo$clusterLabels,
                                      path.expand(outPath))
  invisible(list(path = outPath,
                 variants = as.integer(res$variants),
                 clusters = res$clusters,
                 clusterLabels = clusteringInfo$clusterLabels,
                 clusterSizes = res$clusterSizes))
}

# Resolves the rows of a control matrix file that hold the given case variants.
controlsMatrixVariantRows <- function(info, caseVariants) {
  if (is.null(caseVariants)) {
    return(seq_len(info$nrow))
  }
  rows <- match(as.character(caseVariants), info$rownames)
  if (anyNA(rows)) {
    stop("caseVariants must be a subset of the control matrix variants")
  }
  rows
}
