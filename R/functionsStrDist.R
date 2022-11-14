#' Extract all k-mers of the given string
#'
#' @param str The string to extract k-mers from
#' @param k Size of each k-mer
#'
#' @export
#'
#' @return A string vector containing all k-mers
allKMers <- function(str, k) {
  stopifnot(is.character(str), length(str) <= 1)
  stopifnot(k <= stringr::str_length(str))
  kmers <- c()
  for (i in 1:(stringr::str_length(str)-k+1)) {
    kmers <- c(kmers, substr(str, i, i+k-1))
  }
  return(kmers)
}

countVectors <- function(str1, str2, k) {
  stopifnot(k >= 1)
  stopifnot(is.character(str1), length(str1) <= 1)
  stopifnot(is.character(str2), length(str2) <= 1)
  str1Kmers <- allKMers(str1, k)
  str2Kmers <- allKMers(str2, k)
  levs <- union(str1Kmers, str2Kmers)
  str1Count <- table(factor(str1Kmers, levs))
  str2Count <- table(factor(str2Kmers, levs))
  return(list("first" = str1Count, "second" = str2Count))
}

#' Calculate the Euclidean distance between two
#' strings
#'
#' @param str1 the first string
#' @param str2 the second string
#' @param k size of k-mers
#'
#' @export
#'
#' @return the Euclidean distance between str1 and str2
euclideanDistance <- function(str1, str2, k) {
  vecs <- alfa::countVectors(str1, str2, k)
  stopifnot(length(vecs) == 2)
  str1Count <- vecs$first
  str2Count <- vecs$second
  distMat <- dist(rbind(str1Count, str2Count))
  return(distMat[1])
}

#' Calculate the Euclidean distance between two strings
#' using pre-computed values
#'
fastEuclideanDistance <- function(countVec1, countVec2) {
  return(dist(rbind(countVec1, countVec2))[1])
}

#' Calculate the Euclidean distance between two
#' strings normalized by dividing by the standard deviation
#'
#' @param str1 the first string
#' @param str2 the second string
#' @param k size of k-mers
#'
#' @export
#'
#' @return the standardized Euclidean distance between str1 and str2
standardizedEuclidean <- function(str1, str2, k) {
  vecs <- alfa::standardizedFreq(str1, str2, k)
  stopifnot(length(vecs) == 2)
  str1Freq <- vecs$first
  str2Freq <- vecs$second
  distMat <- dist(rbind(str1Freq, str2Freq))
  return(distMat[1])
}

#' Calculate the standardized Euclidean distance between two
#' strings normalized by dividing by the standard deviation
#'
#' @param str1 the first string
#' @param str2 the second string
#' @param k size of k-mers
#'
#' @export
#'
#' @return the standardized Euclidean distance between str1 and str2
fastStdEuclideanDistance <- function(countVec1, countVec2, str1Len, str2Len,
                                     words, overlapCap, k) {
  freqs1 <- countVec1 / sum(countVec1)
  freqs2 <- countVec2 / sum(countVec2)
  # var1 <- vector(length = length(freqs1))
  # var2 <- vector(length = length(freqs2))
  var1 <- lapply(overlapCap, alfa::freqVariance, seqLen = str1Len, wordLen = k)
  var2 <- lapply(overlapCap, alfa::freqVariance, seqLen = str2Len, wordLen = k)
  dst <- dist(rbind(freqs1 / unlist(var1), freqs2 / unlist(var2)))[1]
  return(dst)
}

normCompressionDistance <- function(str1, str2, k) {

}

#' Calculate the Hamming distance of the two strings
#' (i.e. mean(count vector of str1 != count vector of str2))
#'
#' @param str1 the first string
#' @param str2 the second string
#' @param k size of k-mers
#'
#' @export
#'
#' @return the Hamming distance between str1 and str2
hammingDistance <- function(str1, str2, k) {
  vecs <- alfa::countVectors(str1, str2, k)
  stopifnot(length(vecs) == 2)
  str1Count <- vecs$first
  str2Count <- vecs$second
  return(mean(str1Count != str2Count))
}