readFASTA <- function(filepath) {
  s <- Biostrings::readDNAStringSet(filepath)
  return(s)
}

createDistanceMatrix <- function(dnaSet, metric = "euclidean", k = 10) {
  # stop if input is not of the correct type
  stopifnot(class(dnaSet) == "DNAStringSet")
  # find all combinations of the seq indices
  dnaSetSize <- length(dnaSet)
  pairs <- t(combn(dnaSetSize, 2))
  distMat <- matrix(nrow = dnaSetSize, ncol = dnaSetSize)
  colnames(distMat) <- names(dnaSet)
  rownames(distMat) <- names(dnaSet)
  seqLen <- max(stringr::str_length(dnaSet))

  kmers <- list()
  levs <- c()
  for (i in 1:dnaSetSize) {
    kmer <- alfa::allKMers(as.character(dnaSet[i]), k)
    levs <- union(levs, kmer)
    kmers[[i]] <- kmer
  }
  kmerCounts <- matrix(nrow = dnaSetSize, ncol = length(levs))
  colnames(kmerCounts) <- levs
  for (i in 1:dnaSetSize) {
    kmerCounts[i,] <- table(factor(kmers[[i]], levs))
  }
  if (metric == "standardizedEuclidean")
    overlapCap <- lapply(levs, alfa::overlapCapability)
  for (r in 1:nrow(pairs)) {
    idx1 <- pairs[r, 1]
    idx2 <- pairs[r, 2]
    if (metric == "euclidean")
      dist <- alfa::fastEuclideanDistance(kmerCounts[idx1,], kmerCounts[idx2,])
    else if (metric == "standardizedEuclidean") {
      dist <- alfa::fastStdEuclideanDistance(kmerCounts[idx1,],
                                             kmerCounts[idx2,],
                                             stringr::str_length(dnaSet[idx1]),
                                             stringr::str_length(dnaSet[idx2]),
                                             levs,
                                             overlapCap,
                                             k)
    }

    else
      return(NA)
    stopifnot(!is.na(dist))
    distMat[idx1, idx2] <- dist
    distMat[idx2, idx1] <- dist
  }
  return(distMat)
}

plotPairwiseDist <- function(dnaSet, metric = "euclidean", k = 10, ...) {
  distMat <- alfa::createDistanceMatrix(dnaSet, metric, k)
  plot(distMat, ...)
}
