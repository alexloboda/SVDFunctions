library(testthat)
library(SVDFunctions)

test_that("mvn_stats produces consistent results", {
  # Define covariance matrix, mean vector, test matrix, and clustering in R
  covariance_matrix <- matrix(c(1, 0.8, 0.5,
                                0.8, 1, 0.3,
                                0.5, 0.3, 1), nrow = 3, byrow = TRUE)
  mean_vector <- c(0.9, 0.9, 0.9)
  
  test_matrix <- matrix(c(1.46611034, 2.1128962, 0.6487757,
                          1.18410660, 1.3453171, 1.0116327,
                         -0.71017476, 0.3972155, -0.6873396,
                          1.69569546, 0.6755756, 0.2220538,
                          1.33330116, 0.8317992, 0.3448378,
                         -0.72172455, 0.4980769, -1.2921344,
                          0.32584272, 0.3555474, 1.3457023,
                          2.02847892, 2.1606610, 2.0010079,
                          2.04383963, 1.8677370, 0.5502761,
                          1.21905946, 0.8897409, 2.2307804,
                         -0.13594810, -0.3100018, 0.4732567,
                          0.69358431, 0.8163033, 0.5440840,
                          0.48405706, 0.2903019, 1.3962302,
                          1.18934191, 0.1821820, 1.4522213,
                          1.34002737, 1.2023193, 1.9982008,
                         -0.31851695, -1.2387499, 0.1778412,
                         -0.01797021, 0.7399531, 1.2134355,
                          3.03958306, 2.5874102, 2.2179830,
                          0.41968262, 0.4833934, 0.2940008,
                          1.56444597, 1.5137194, 1.0163329), nrow = 20, byrow = TRUE)

  clustering <- c(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6)
  # Call the combined C++ function to get stats
  stats <- SVDFunctions:::rcpp_run_mvn_stats_tests_combined(covariance_matrix, mean_vector, test_matrix, clustering)
  
  transform <- function(l) {
    m <- matrix(unlist(l), nrow = 7)
  }
  
  stats$interpoint_pairwise <- transform(stats$interpoint_pairwise)
  stats$approx_pairwise <- transform(stats$approx_pairwise)
  
  # Perform checks
  centered_stats_diff <- abs(stats$interpoint_centered - stats$approx_centered)
  expect_true(all(centered_stats_diff <= 0.05), "Centered stats comparison failed")
  
  pairwise_stats_diff <- abs(stats$interpoint_pairwise - stats$approx_pairwise)
  expect_true(all(pairwise_stats_diff <= 0.05), "Pairwise stats comparison failed")
  #print maximum pairwise difference
  max_pairwise_diff <- max(pairwise_stats_diff)
  print(paste("Maximum pairwise stats difference:", max_pairwise_diff))
})

test_that("mvn_stats works with 5D, 5 clusters, synthetic Gaussian data", {
  set.seed(42)
  library(MASS)
  d <- 5
  n_per_cluster <- 25
  k <- 5
  n <- n_per_cluster * k

  # Covariance matrix: some off-diagonal zeros, but not all
  covariance_matrix <- matrix(c(
    1,   0.5, 0,   0,   0.2,
    0.5, 1,   0.3, 0,   0,
    0,   0.3, 1,   0.4, 0,
    0,   0,   0.4, 1,   0.1,
    0.2, 0,   0,   0.1, 1
  ), nrow = d, byrow = TRUE)

  # Cluster means
  mean_vectors <- matrix(runif(d * k, min = -2, max = 2), nrow = k)

  # Generate data
  test_matrix <- do.call(rbind, lapply(1:k, function(i) {
    MASS::mvrnorm(n_per_cluster, mu = mean_vectors[i,], Sigma = covariance_matrix)
  }))

  clustering <- rep(0:(k-1), each = n_per_cluster)

  stats <- SVDFunctions:::rcpp_run_mvn_stats_tests_combined(covariance_matrix, colMeans(test_matrix), test_matrix, clustering)

  transform <- function(l) {
    matrix(unlist(l), nrow = k)
  }
  stats$interpoint_pairwise <- transform(stats$interpoint_pairwise)
  stats$approx_pairwise <- transform(stats$approx_pairwise)

  centered_stats_diff <- abs(stats$interpoint_centered - stats$approx_centered)
  #expect_true(all(centered_stats_diff <= 0.05), "Centered stats comparison failed (5D synthetic)")

  pairwise_stats_diff <- abs(stats$interpoint_pairwise - stats$approx_pairwise)
  #expect_true(all(pairwise_stats_diff <= 0.05), "Pairwise stats comparison failed (5D synthetic)")
  
  max_abs_stat <- max(abs(stats$interpoint_pairwise))
  max_rel_error <- pairwise_stats_diff / max_abs_stat

  max_pairwise_diff <- max(pairwise_stats_diff)
  print(paste("Maximum pairwise stats difference (5D synthetic):", max_pairwise_diff))
  print(paste("Maximum relative error (5D synthetic):", max(max_rel_error)))
})
