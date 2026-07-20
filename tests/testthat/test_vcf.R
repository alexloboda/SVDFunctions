context("parsing VCF files")

vcf <- "CEU.exon.2010_09.genotypes.vcf.gz"
file <- system.file("extdata", vcf, package = "SVDFunctions")
DP <- 20

missing_function <- function(x) {
  missing_function_dp(x, DP)
}

missing_function_dp <- function(x, dp) {
  vapply(x, function(value) {
    is.na(value) || startsWith(value, ".") ||
      as.integer(strsplit(value, ":", fixed = TRUE)[[1]][2]) < dp
  }, logical(1))
}

allele_string <- function(x) {
  if (is.logical(x)) {
    ifelse(is.na(x), NA_character_, ifelse(x, "T", "F"))
  } else {
    as.character(x)
  }
}

parse_gt_value <- function(x, dp) {
  if (missing_function_dp(x, dp)) {
    return(NA_real_)
  }
  if (startsWith(x, "0/1") || startsWith(x, "1/0")) {
    return(1)
  }
  if (startsWith(x, "0/0")) {
    return(0)
  }
  if (startsWith(x, "1/1")) {
    return(2)
  }
  NA_real_
}

extract_genotype_matrix <- function(df, samples, dp) {
  gt <- apply(df[, samples, drop = FALSE], c(1, 2), parse_gt_value, dp = dp)
  if (!is.matrix(gt)) {
    gt <- matrix(gt, nrow = 1)
  }
  colnames(gt) <- samples
  gt
}

# Mirrors src/vcf_r.cpp::Counts::pass to keep the direct VCF ground truth
# aligned with the binary scanner's filtering semantics.
counts_pass_binary_filters <- function(counts, n_samples,
                                       minMAF = 0, maxMAF = 1,
                                       minCallRate = 0.9,
                                       minMAC = 1L,
                                       maxMAC = .Machine$integer.max,
                                       reportSingletons = TRUE) {
  eps <- 1e-6
  left <- 2 * counts[1] + counts[2]
  right <- 2 * counts[3] + counts[2]
  sum_counts <- sum(counts)
  if (left < right) {
    tmp <- left
    left <- right
    right <- tmp
  }
  maf <- right / (left + right)
  call_rate <- sum_counts / n_samples
  singletons <- reportSingletons || (counts[1] + counts[2] != 1 && counts[2] + counts[3] != 1)
  call_rate > minCallRate &&
    maf + eps > minMAF &&
    maf - eps < maxMAF &&
    right >= minMAC &&
    right <= maxMAC &&
    left > 0 &&
    right > 0 &&
    singletons
}

direct_variant_counts_from_vcf <- function(vcfFile, request, samples, dp,
                                           minMAF = 0, maxMAF = 1,
                                           minCallRate = 0.9,
                                           minMAC = 1L,
                                           maxMAC = .Machine$integer.max,
                                           reportSingletons = TRUE) {
  tokens <- strsplit(request, "\t", fixed = TRUE)[[1]]
  if (length(tokens) != 3L) {
    return(NULL)
  }
  position_tokens <- strsplit(tokens[1], ":", fixed = TRUE)[[1]]
  if (length(position_tokens) != 2L) {
    return(NULL)
  }
  query <- paste0(sub("^chr", "", position_tokens[1]), ":", position_tokens[2], "-", position_tokens[2])
  df <- seqminer::tabix.read.table(vcfFile, query)
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }
  labels <- paste0("chr", df$CHROM, ":", df$POS, "\t", allele_string(df$REF), "\t", allele_string(df$ALT))
  hit <- match(request, labels)
  reversed <- FALSE
  if (is.na(hit)) {
    reversed_request <- paste(tokens[1], tokens[3], tokens[2], sep = "\t")
    hit <- match(reversed_request, labels)
    reversed <- !is.na(hit)
  }
  if (is.na(hit)) {
    return(NULL)
  }

  gt <- extract_genotype_matrix(df[hit, , drop = FALSE], samples, dp)
  counts <- genotypesToCounts(gt)[1, ]
  if (!counts_pass_binary_filters(counts, length(samples), minMAF, maxMAF,
                                  minCallRate, minMAC, maxMAC, reportSingletons)) {
    return(NULL)
  }
  if (reversed) {
    counts <- counts[c("hom_alt", "het", "hom_ref")]
    names(counts) <- c("hom_ref", "het", "hom_alt")
  }
  matrix(counts, nrow = 1, dimnames = list(request, c("hom_ref", "het", "hom_alt")))
}

direct_region_counts_from_vcf <- function(vcfFile, region, samples, dp,
                                          minMAF = 0, maxMAF = 1,
                                          minCallRate = 0.9,
                                          minMAC = 1L,
                                          maxMAC = .Machine$integer.max,
                                          reportSingletons = TRUE) {
  tokens <- strsplit(region, " +")[[1]]
  query <- paste0(sub("^chr", "", tokens[1]), ":", tokens[2], "-", tokens[3])
  df <- seqminer::tabix.read.table(vcfFile, query)
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }

  gt <- extract_genotype_matrix(df, samples, dp)
  counts <- genotypesToCounts(gt)
  keep <- apply(
    counts,
    1,
    counts_pass_binary_filters,
    n_samples = length(samples),
    minMAF = minMAF,
    maxMAF = maxMAF,
    minCallRate = minCallRate,
    minMAC = minMAC,
    maxMAC = maxMAC,
    reportSingletons = reportSingletons
  )
  if (!any(keep)) {
    return(NULL)
  }
  matrix(colSums(counts[keep, , drop = FALSE]),
         nrow = 1,
         dimnames = list(region, c("hom_ref", "het", "hom_alt")))
}

bind_non_null_matrices <- function(items) {
  items <- Filter(Negate(is.null), items)
  if (!length(items)) {
    matrix(numeric(0), nrow = 0, ncol = 3,
           dimnames = list(NULL, c("hom_ref", "het", "hom_alt")))
  } else {
    do.call(rbind, items)
  }
}

test_that("genotype matrix are parsed correctly", {
  bannedPos <- "chr1:18022097"
  banned <- "chr1:18022097\tG\tT"
  samples <- sampleNamesVCF(file)
  vcf <- genotypeMatrixVCF(file, DP = DP, GQ = 0, samples = samples[1:10], 
                               bannedPositions = bannedPos)
  df <- seqminer::tabix.read.table(file, paste0(c(1:22), ":1-300000000"))
  format <- function(x) paste0("chr", df$CHROM[x], ":", df$POS[x], "\t", 
                               df$REF[x], "\t", df$ALT[x])
  rownames(df) <- lapply(1:nrow(df), format)
  df <- df[, samples[1:10]]
  df <- apply(df, c(1, 2), function(x) if(missing_function(x)) NA else x)
  df <- apply(df, c(1, 2), function(x) {
    if (is.na(x)) {
      NA
    } else if (startsWith(x, "0/1") | startsWith(x, "1/0")) {
      1
    } else if (startsWith(x, "0/0")) {
      0
    } else if (startsWith(x, "1/1")) {
      2
    }
  })
  df <- df[rowSums(is.na(df)) <= 1, ]
  df <- df[rowSums(matrix(!is.na(df) & df == 1, nrow = nrow(df))) > 0 | 
          (rowSums(matrix(!is.na(df) & df == 0, nrow = nrow(df))) > 0 & 
           rowSums(matrix(!is.na(df) & df == 2, nrow = nrow(df))) > 0), ]
  df <- df[rownames(df) != banned, ]
  expect_equal(df, vcf)
})

test_that("inversed variants are read correctly", {
  var <- "chr1:18022153\tA\tC"
  samples <- sampleNamesVCF(file)[1:10]
  vcf <- genotypeMatrixVCF(file, DP = DP, GQ = 0, samples = samples, 
                           variants = var)
  expected <- matrix(c(0, 1, 0, 2, 1, 2, NaN, 1, 0, 0),
                     nrow = 1) 
  rownames(expected) <- var
  colnames(expected) <- samples
  expect_equal(expected, vcf)
})

test_that("callrates are calculated correctly", {
  regions <- data.frame(chr = c("1", "1", "20"), 
                        from = c("1108138", "40000000", "33521213"), 
                        to = c("36000000", "40000010", "61665700"))
  samples <- c("NA06989", "NA10847", "NA11840", "NA12873")
  pkgFormat = function(x) paste0("chr", x[1], " ", x[2], " ", x[3])
  seqMinerFormat = function(x) paste0(x[1], ":", x[2], "-", x[3])
  regionsPkg <- apply(regions, 1, pkgFormat)
  regionsSeqMiner <- apply(regions, 1, seqMinerFormat)
  cr <- callRateMatrixVCF(file, DP = DP, GQ = 0, regions = regionsPkg, 
                           samples = samples)
  expectedMatrix <- matrix(nrow = 0, ncol = length(samples))
  for (i in 1:nrow(regions)) {
    df <- seqminer::tabix.read.table(file, regionsSeqMiner[i])
    if (nrow(df) == 0) {
      expectedMatrix <- rbind(expectedMatrix, NA)
      next
    }
    df <- df[, samples]
    missing <- sapply(df[, samples], function(x) sum(missing_function(x)))
    expected <- (nrow(df) - missing) / nrow(df)
    expectedMatrix <- rbind(expectedMatrix, expected)
  }
  rownames(expectedMatrix) <- regionsPkg
  expectedMatrix <- expectedMatrix[apply(expectedMatrix, 1, 
                                         function(x) !any(is.na(x))), ]
  expect_equal(cr, expectedMatrix, tolerance = 1e-8)
})

test_that("storing/extracting data to/from binary file works ", {
  set.seed(42)
  prefix <- paste0(tempdir(), "/db")
  bin <- paste0(prefix, "_bin")
  meta <- paste0(prefix, "_meta")
  vcf <- scanVCF(file, DP = DP, GQ = 0, binaryPath = prefix)
  
  samples <- vcf$samples
  variants <- rownames(vcf$genotype)
  variants <- c(variants, "chr1:1\nT\nC")
  samples <- sample(samples, as.integer(length(samples) / 2))
  sample_size <- length(variants) %/% 2
  sel <- c(rep(TRUE, sample_size), rep(FALSE, length(variants) - sample_size))
  variants <- variants[sample(sel)]
  
  selReverse <- sample(length(variants) %/% 3)
  variants[selReverse] <- sapply(variants[selReverse], function(x) {
    tokens <- unlist(strsplit(x, "\t"))
    paste(tokens[1], tokens[3], tokens[2], sep = "\t")
  })
  
  localDP <- 30
  actual <- scanBinaryFile(bin, meta, samples, variants, DP = localDP, GQ = 0)
  expected <- bind_non_null_matrices(lapply(
    variants,
    direct_variant_counts_from_vcf,
    vcfFile = file,
    samples = samples,
    dp = localDP
  ))
  expect_equal(actual[, 1:3], expected)
  
  regions <- c("chr2 148943492 234267016",
               "chr7 102949810 149848267",
               "chr10 46506985 46507378")
  
  expectedReg <- list()
  expectedReg[["complete"]] <- bind_non_null_matrices(lapply(
    regions,
    direct_region_counts_from_vcf,
    vcfFile = file,
    samples = samples,
    dp = localDP
  ))
  expectedReg[["maf"]] <- bind_non_null_matrices(lapply(
    regions,
    direct_region_counts_from_vcf,
    vcfFile = file,
    samples = samples,
    dp = localDP,
    minMAF = 0.04
  ))
  expectedReg[["rare"]] <- bind_non_null_matrices(lapply(
    regions,
    direct_region_counts_from_vcf,
    vcfFile = file,
    samples = samples,
    dp = localDP,
    maxMAF = 0.04
  ))
  expectedReg[["cr"]] <- bind_non_null_matrices(lapply(
    regions,
    direct_region_counts_from_vcf,
    vcfFile = file,
    samples = samples,
    dp = localDP,
    minCallRate = 41.5 / 45
  ))
  
  reg <- list()
  reg[["complete"]] <- scanBinaryFile(bin, meta, samples, regions = regions, 
                        DP = localDP, GQ = 0)[, 1:3, drop = FALSE]
  reg[["maf"]] <- scanBinaryFile(bin, meta, samples, regions = regions, 
                        DP = localDP, GQ = 0, minMAF = 0.04)[, 1:3, drop = FALSE]
  reg[["rare"]] <- scanBinaryFile(bin, meta, samples, regions = regions, 
                        DP = localDP, GQ = 0, maxMAF = 0.04)[, 1:3, drop = FALSE]
  reg[["cr"]] <- scanBinaryFile(bin, meta, samples, regions = regions, 
                        DP = localDP, GQ = 0, minCallRate = 41.5 / 45)[, 1:3, drop = FALSE]
  
  expect_equal(expectedReg, reg)
})

test_that("parsing multivariant lines works", {
  file <- system.file("extdata", "multivariant.vcf.gz", package = "SVDFunctions")
  GT <- genotypeMatrixVCF(file, DP = 0, GQ = 0)
  expected <- matrix(c(0, 1, 1, 1, 1, 2), ncol = 3)
  colnames(expected) <- c("A", "B", "C")
  rownames(expected) <- c("chr1:1\tT\tG", "chr1:2\tT\t*")
  expect_equal(GT, expected)
})

# Helper that splits sample names into k clusters in a round-robin fashion and
# returns a named vector (names = samples, values = cluster labels).
assign_round_robin_clusters <- function(samples, k) {
  labels <- paste0("cl", seq_len(k))
  assignment <- labels[((seq_along(samples) - 1) %% k) + 1]
  names(assignment) <- samples
  assignment
}

test_that("per-cluster binary aggregation matches sample-level scan", {
  set.seed(7)
  srcPrefix <- paste0(tempdir(), "/cluster_src")
  srcBin <- paste0(srcPrefix, "_bin")
  srcMeta <- paste0(srcPrefix, "_meta")

  # Canonical source binary keeps every genotype (DP = 0, GQ = 0) so that the
  # DP/GQ thresholds applied at conversion time are exactly reproducible.
  vcf <- scanVCF(file, DP = 0, GQ = 0, binaryPathPrefix = srcPrefix)
  samples <- vcf$samples

  clusters <- assign_round_robin_clusters(samples, 4)
  clusterPrefix <- paste0(tempdir(), "/cluster_dst")
  buildDP <- 20L

  info <- buildClusterBinaryFromBinary(srcBin, srcMeta, clusters, clusterPrefix,
                                       DP = buildDP, GQ = 0L)
  clBin <- info$binaryFile
  clMeta <- info$metafile

  expect_true(file.exists(clBin))
  expect_true(file.exists(clMeta))
  expect_equal(sort(info$clusters), sort(unique(clusters)))
  # On-disk layout: variants x clusters x sizeof(ClusterCounts) (3 bytes).
  expect_equal(file.size(clBin),
               info$variants * length(info$clusters) * 3)

  regions <- c("chr2 148943492 234267016",
               "chr7 102949810 149848267",
               "chr10 46506985 46507378")

  selectedClusters <- c("cl1", "cl3")
  selectedSamples <- names(clusters)[clusters %in% selectedClusters]

  for (mcr in c(0, 0.9)) {
    actual <- scanClusterBinaryFile(clBin, clMeta, clusters = selectedClusters,
                                    regions = regions, minCallRate = mcr)
    expected <- scanBinaryFile(srcBin, srcMeta, selectedSamples,
                               regions = regions, DP = buildDP, GQ = 0,
                               minCallRate = mcr)
    expect_equal(actual, expected)
  }
})

test_that("per-cluster scan reproduces variant-level counts and MAF filters", {
  srcPrefix <- paste0(tempdir(), "/cluster_src2")
  srcBin <- paste0(srcPrefix, "_bin")
  srcMeta <- paste0(srcPrefix, "_meta")

  vcf <- scanVCF(file, DP = 0, GQ = 0, binaryPathPrefix = srcPrefix)
  samples <- vcf$samples
  storedVariants <- rownames(vcf$genotype)

  clusters <- assign_round_robin_clusters(samples, 3)
  clusterPrefix <- paste0(tempdir(), "/cluster_dst2")
  buildDP <- 20L
  info <- buildClusterBinaryFromBinary(srcBin, srcMeta, clusters, clusterPrefix,
                                       DP = buildDP, GQ = 0L)

  reqVariants <- head(storedVariants, 25)
  selectedClusters <- c("cl2", "cl3")
  selectedSamples <- names(clusters)[clusters %in% selectedClusters]

  actual <- scanClusterBinaryFile(info$binaryFile, info$metafile,
                                  clusters = selectedClusters,
                                  variants = reqVariants, minCallRate = 0,
                                  minMAF = 0.04)
  expected <- scanBinaryFile(srcBin, srcMeta, selectedSamples,
                             variants = reqVariants, DP = buildDP, GQ = 0,
                             minCallRate = 0, minMAF = 0.04)
  expect_equal(actual, expected)
})

test_that("per-cluster scan with all clusters equals full sample scan", {
  srcPrefix <- paste0(tempdir(), "/cluster_src3")
  srcBin <- paste0(srcPrefix, "_bin")
  srcMeta <- paste0(srcPrefix, "_meta")

  vcf <- scanVCF(file, DP = 0, GQ = 0, binaryPathPrefix = srcPrefix)
  samples <- vcf$samples

  clusters <- assign_round_robin_clusters(samples, 5)
  clusterPrefix <- paste0(tempdir(), "/cluster_dst3")
  buildDP <- 20L
  info <- buildClusterBinaryFromBinary(srcBin, srcMeta, clusters, clusterPrefix,
                                       DP = buildDP, GQ = 0L)

  regions <- c("chr2 148943492 234267016", "chr7 102949810 149848267")

  # clusters = NULL aggregates every cluster, i.e. every sample.
  actual <- scanClusterBinaryFile(info$binaryFile, info$metafile,
                                  clusters = NULL, regions = regions,
                                  minCallRate = 0)
  expected <- scanBinaryFile(srcBin, srcMeta, samples, regions = regions,
                             DP = buildDP, GQ = 0, minCallRate = 0)
  expect_equal(actual, expected)
})

test_that("per-cluster binary input is validated", {
  srcPrefix <- paste0(tempdir(), "/cluster_src4")
  srcBin <- paste0(srcPrefix, "_bin")
  srcMeta <- paste0(srcPrefix, "_meta")

  vcf <- scanVCF(file, DP = 0, GQ = 0, binaryPathPrefix = srcPrefix)
  samples <- vcf$samples
  clusters <- assign_round_robin_clusters(samples, 3)
  clusterPrefix <- paste0(tempdir(), "/cluster_dst4")
  info <- buildClusterBinaryFromBinary(srcBin, srcMeta, clusters, clusterPrefix,
                                       DP = 20L, GQ = 0L)

  # Unnamed clustering is rejected.
  expect_error(
    buildClusterBinaryFromBinary(srcBin, srcMeta, unname(clusters),
                                 paste0(tempdir(), "/cluster_dst_bad"),
                                 DP = 20L, GQ = 0L),
    "named vector"
  )

  # Sample absent from the metadata is rejected.
  badClusters <- clusters
  names(badClusters)[1] <- "definitely_not_a_sample"
  expect_error(
    buildClusterBinaryFromBinary(srcBin, srcMeta, badClusters,
                                 paste0(tempdir(), "/cluster_dst_bad2"),
                                 DP = 20L, GQ = 0L),
    "not found"
  )

  # Unknown cluster requested at scan time is rejected.
  expect_error(
    scanClusterBinaryFile(info$binaryFile, info$metafile,
                          clusters = "no_such_cluster",
                          regions = "chr2 148943492 234267016"),
    "not found"
  )

  # Truncated cluster-size metadata is rejected.
  brokenMeta <- paste0(tempdir(), "/cluster_dst4_broken_meta")
  brokenLines <- readLines(info$metafile)
  brokenLines[2] <- sub("\t.*$", "\t", brokenLines[2])
  writeLines(brokenLines, brokenMeta)
  expect_error(
    scanClusterBinaryFile(info$binaryFile, brokenMeta,
                          regions = "chr2 148943492 234267016",
                          minCallRate = 0),
    "incomplete cluster sizes"
  )
})

test_that("buildClusterBinaryFromBinary warns about dropped samples", {
  srcPrefix <- paste0(tempdir(), "/cluster_src5")
  srcBin <- paste0(srcPrefix, "_bin")
  srcMeta <- paste0(srcPrefix, "_meta")

  vcf <- scanVCF(file, DP = 0, GQ = 0, binaryPathPrefix = srcPrefix)
  samples <- vcf$samples

  # Covering only a subset of the source samples emits a warning.
  subset <- samples[seq_len(length(samples) %/% 2)]
  partialClusters <- assign_round_robin_clusters(subset, 2)
  expect_warning(
    buildClusterBinaryFromBinary(srcBin, srcMeta, partialClusters,
                                 paste0(tempdir(), "/cluster_dst5"),
                                 DP = 20L, GQ = 0L),
    "dropped"
  )

  # Full coverage does not warn.
  fullClusters <- assign_round_robin_clusters(samples, 2)
  expect_no_warning(
    buildClusterBinaryFromBinary(srcBin, srcMeta, fullClusters,
                                 paste0(tempdir(), "/cluster_dst5b"),
                                 DP = 20L, GQ = 0L)
  )
})

test_that("indels are matched as-is while SNVs may be flipped (binary + cluster)", {
  # The shared CEU test VCF contains only SNVs, so a bespoke biallelic VCF is
  # needed to exercise indel handling. It holds one SNV, one insertion and one
  # deletion, each with both alleles observed (so MAC_filter keeps them) and
  # counts hom_ref = 2, het = 1, hom_alt = 1.
  samples4 <- c("S1", "S2", "S3", "S4")
  vcfLines <- c(
    "##fileformat=VCFv4.1",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", samples4), collapse = "\t"),
    paste("1", "100", ".", "A", "G", ".", "PASS", ".", "GT:DP",
          "0/0:30", "0/1:30", "1/1:30", "0/0:30", sep = "\t"),
    paste("1", "200", ".", "C", "CAT", ".", "PASS", ".", "GT:DP",
          "0/1:30", "0/0:30", "0/0:30", "1/1:30", sep = "\t"),
    paste("1", "300", ".", "GTT", "G", ".", "PASS", ".", "GT:DP",
          "1/1:30", "0/1:30", "0/0:30", "0/0:30", sep = "\t")
  )
  vcfPath <- paste0(tempdir(), "/indel_synth.vcf.gz")
  con <- gzfile(vcfPath, "w")
  writeLines(vcfLines, con)
  close(con)

  prefix <- paste0(tempdir(), "/indel_synth")
  srcBin <- paste0(prefix, "_bin")
  srcMeta <- paste0(prefix, "_meta")
  scanVCF(vcfPath, DP = 0, GQ = 0, binaryPathPrefix = prefix)

  snvAsIs <- "chr1:100\tA\tG"
  snvFlip <- "chr1:100\tG\tA"
  insAsIs <- "chr1:200\tC\tCAT"
  insFlip <- "chr1:200\tCAT\tC"
  delAsIs <- "chr1:300\tGTT\tG"
  delFlip <- "chr1:300\tG\tGTT"

  scanBin <- function(variants) {
    scanBinaryFile(srcBin, srcMeta, samples4, variants = variants,
                   DP = 0, GQ = 0, minCallRate = 0, minMAC = 0L)
  }

  # A flipped SNV is answered with the stored variant, reported under the
  # requested (reversed) name and with hom_ref/hom_alt swapped.
  asIs <- scanBin(snvAsIs)
  flip <- scanBin(snvFlip)
  expect_equal(rownames(asIs), snvAsIs)
  expect_equal(rownames(flip), snvFlip)
  expect_equal(unname(asIs[, c("hom_ref", "het", "hom_alt")]), c(2, 1, 1))
  expect_equal(unname(flip[, "hom_ref"]), unname(asIs[, "hom_alt"]))
  expect_equal(unname(flip[, "hom_alt"]), unname(asIs[, "hom_ref"]))
  expect_equal(unname(flip[, "het"]), unname(asIs[, "het"]))

  # Indels are matched in the stored orientation only: swapping their alleles
  # describes a different event, so a reversed request must not match.
  expect_equal(rownames(scanBin(insAsIs)), insAsIs)
  expect_equal(rownames(scanBin(delAsIs)), delAsIs)
  expect_equal(nrow(scanBin(insFlip)), 0L)
  expect_equal(nrow(scanBin(delFlip)), 0L)

  # The per-cluster scanner must show the same semantics (flipped SNVs present,
  # flipped indels absent).
  clusters4 <- c(S1 = "clA", S2 = "clA", S3 = "clB", S4 = "clB")
  clPrefix <- paste0(tempdir(), "/indel_synth_cl")
  info <- buildClusterBinaryFromBinary(srcBin, srcMeta, clusters4, clPrefix,
                                       DP = 0L, GQ = 0L)

  scanCl <- function(variants) {
    scanClusterBinaryFile(info$binaryFile, info$metafile,
                          clusters = c("clA", "clB"), variants = variants,
                          minCallRate = 0, minMAC = 0L)
  }

  clAsIs <- scanCl(snvAsIs)
  clFlip <- scanCl(snvFlip)
  expect_equal(rownames(clAsIs), snvAsIs)
  expect_equal(rownames(clFlip), snvFlip)
  expect_equal(unname(clFlip[, "hom_ref"]), unname(clAsIs[, "hom_alt"]))
  expect_equal(unname(clFlip[, "hom_alt"]), unname(clAsIs[, "hom_ref"]))
  expect_equal(unname(clFlip[, "het"]), unname(clAsIs[, "het"]))
  expect_equal(nrow(scanCl(insFlip)), 0L)
  expect_equal(nrow(scanCl(delFlip)), 0L)

  # Binary and cluster scans agree on the flipped SNV counts.
  expect_equal(clFlip[, c("hom_ref", "het", "hom_alt")],
               flip[, c("hom_ref", "het", "hom_alt")])
})

