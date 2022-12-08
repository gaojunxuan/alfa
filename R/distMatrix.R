#' Read FASTA file into a DNAStringSet object to be processed by other
#' functions.
#'
#' @param filepath path to the FASTA file to be imported
#'
#' @examples
#' alfa::readFASTA("inst/extdata/usflu.fasta")
#'
#' @return a DNAStringSet object
#'
#' @export
readFASTA <- function(filepath) {
  s <- Biostrings::readDNAStringSet(filepath)
  return(s)
}

#' Create a pairwise distance matrix for the given set of strings
#'
#' Although the parameter is named dnaSet, it can be any set of strings of
#' the type DNAStringSet and AAStringSet over any well-defined alphabet
#'
#' @param dnaSet an object of type DNAStringSet or AAStringSet
#' @param metric metric for computing the distance, for now only 'euclidean'
#' and 'standardizedEuclidean' are supported
#' @param k size of k-mers (10 by default for DNA strings)
#'
#' @return an length(dnaSet) by length(dnaSet) matrix containing the pairwise
#' distance of every pair of strings in the set with the diagonal entries set
#' to NA.
#'
#' @examples
#' library("ape")
#' data("woodmouse")
#' woodmouse <- as.list(woodmouse)
#' woodmouse <- lapply(woodmouse, ape::as.character.DNAbin)
#' woodmouse <- lapply(woodmouse, paste0, collapse = "")
#' woodmouseSet <- Biostrings::DNAStringSet(unlist(woodmouse))
#' alfa::createDistanceMatrix(woodmouseSet, metric = "euclidean")
#'
#' @import Biostrings
#'
#' @export
createDistanceMatrix <- function(dnaSet, metric = "euclidean", k = 10) {
  # stop if input is not of the correct type
  stopifnot(class(dnaSet) == "DNAStringSet" || class(dnaSet) == "AAStringSet")
  # find all combinations of the seq indices
  dnaSetSize <- length(dnaSet)
  pairs <- t(combn(dnaSetSize, 2))
  distMat <- matrix(nrow = dnaSetSize, ncol = dnaSetSize)
  colnames(distMat) <- names(dnaSet)
  rownames(distMat) <- names(dnaSet)
  seqLen <- max(stringr::str_length(dnaSet))

  # precompute the level vectors and kmer counts
  kmers <- list()
  levs <- c()
  for (i in 1:dnaSetSize) {
    kmer <- allKMers(as.character(dnaSet[i]), k)
    levs <- union(levs, kmer)
    kmers[[i]] <- kmer
  }
  kmerCounts <- matrix(nrow = dnaSetSize, ncol = length(levs))
  colnames(kmerCounts) <- levs
  for (i in 1:dnaSetSize) {
    kmerCounts[i,] <- table(factor(kmers[[i]], levs))
  }
  # compute overlapping capability for stdEuclidean
  if (metric == "standardizedEuclidean") {
    overlapCap <- lapply(levs, overlapCapability)
  }
  for (r in 1:nrow(pairs)) {
    idx1 <- pairs[r, 1]
    idx2 <- pairs[r, 2]
    if (metric == "euclidean") {
      dist <- fastEuclideanDistance(kmerCounts[idx1,], kmerCounts[idx2,])
    }
    else if (metric == "standardizedEuclidean") {
      dist <- fastStdEuclideanDistance(kmerCounts[idx1,],
                                             kmerCounts[idx2,],
                                             stringr::str_length(dnaSet[idx1]),
                                             stringr::str_length(dnaSet[idx2]),
                                             levs,
                                             overlapCap,
                                             k)
    }
    # not a supported metric
    else {
      return(NA)
    }
    stopifnot(!is.na(dist))
    distMat[idx1, idx2] <- dist
    distMat[idx2, idx1] <- dist
  }
  return(distMat)
}

#' Create and then plot the distance matrix as a heatmap
#' @param dnaSet an object of type DNAStringSet or AAStringSet
#' @param metric metric for computing the distance, for now only 'euclidean'
#' and 'standardizedEuclidean' are supported
#' @param k size of k-mers (10 by default for DNA strings)
#' @param ... arguments to be passed to the plotting function 'heatmap'
#'
#' @examples
#' library("ape")
#' data("woodmouse")
#' woodmouse <- as.list(woodmouse)
#' woodmouse <- lapply(woodmouse, ape::as.character.DNAbin)
#' woodmouse <- lapply(woodmouse, paste0, collapse = "")
#' woodmouseSet <- Biostrings::DNAStringSet(unlist(woodmouse))
#' alfa::plotPairwiseDist(woodmouseSet, metric = "euclidean")
#'
#' @import Biostrings
#'
#' @export
plotPairwiseDist <- function(dnaSet, metric = "euclidean", k = 10, ...) {
  distMat <- createDistanceMatrix(dnaSet, metric, k)
  heatmap(distMat, ...)
}
