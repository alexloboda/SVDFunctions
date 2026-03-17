exactWideSvd <- function(x, targetRank, nu, nv) {
  basisRank <- max(nu, nv)
  gram <- tcrossprod(x)
  eigenResult <- eigen(gram, symmetric = TRUE)

  singularValues <- sqrt(pmax(eigenResult$values[seq_len(targetRank)], 0))
  leftBasis <- eigenResult$vectors[, seq_len(basisRank), drop = FALSE]

  rightBasis <- matrix(numeric(0), nrow = ncol(x), ncol = 0L)
  if (nv > 0L) {
    rightBasis <- matrix(0, nrow = ncol(x), ncol = nv)
    stable <- singularValues[seq_len(nv)] > (.Machine$double.eps^0.5) * max(1, singularValues[1])
    if (any(stable)) {
      stableIds <- which(stable)
      projected <- crossprod(x, leftBasis[, stableIds, drop = FALSE])
      rightBasis[, stableIds] <- sweep(projected, 2L, singularValues[stableIds], "/")
    }
  }

  list(
    u = leftBasis[, seq_len(nu), drop = FALSE],
    d = singularValues,
    v = rightBasis
  )
}

truncatedSvd <- function(x, k, nu = k, nv = k, work = NULL) {
  dims <- dim(x)

  if (length(dims) != 2L) {
    stop("Expected a matrix for SVD.", call. = FALSE)
  }

  maxRank <- min(dims)
  if (maxRank < 1L) {
    stop("Cannot compute SVD of an empty matrix.", call. = FALSE)
  }

  k <- as.integer(k)[1]
  nu <- as.integer(nu)[1]
  nv <- as.integer(nv)[1]

  if (is.na(k) || k < 1L) {
    stop("k must be a positive integer.", call. = FALSE)
  }

  if (is.na(nu) || nu < 0L) {
    stop("nu must be a non-negative integer.", call. = FALSE)
  }

  if (is.na(nv) || nv < 0L) {
    stop("nv must be a non-negative integer.", call. = FALSE)
  }

  if (!is.null(work)) {
    work <- as.integer(work)[1]
    if (is.na(work) || work < 1L) {
      stop("work must be a positive integer.", call. = FALSE)
    }
  }

  targetRank <- min(k, maxRank)
  nu <- min(nu, nrow(x), targetRank)
  nv <- min(nv, ncol(x), targetRank)
  useExactWideFallback <- ncol(x) > nrow(x)

  if (targetRank == maxRank || (2L * targetRank) >= maxRank) {
    if (useExactWideFallback) {
      return(exactWideSvd(x, targetRank = targetRank, nu = nu, nv = nv))
    }

    svdResult <- base::svd(x, nu = nu, nv = nv)
    svdResult$d <- svdResult$d[seq_len(targetRank)]
    return(svdResult)
  }

  if (!is.null(work)) {
    work <- max(work, targetRank + 7L)
  }

  svdResult <- tryCatch(
    if (is.null(work)) {
      irlba::irlba(x, nu = targetRank, nv = targetRank)
    } else {
      irlba::irlba(x, nu = targetRank, nv = targetRank, work = work)
    },
    error = function(err) {
      if (useExactWideFallback) {
        return(exactWideSvd(x, targetRank = targetRank, nu = nu, nv = nv))
      }

      stop(err)
    }
  )
  list(
    u = svdResult$u[, seq_len(nu), drop = FALSE],
    d = svdResult$d[seq_len(targetRank)],
    v = svdResult$v[, seq_len(nv), drop = FALSE]
  )
}