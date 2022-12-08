#' Extract all k-mers of the given string
#'
#' @param str The string to extract k-mers from
#' @param k Size of each k-mer
#'
#' @examples
#' alfa::allKMers("ATGTGAC", 3)
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

#' Create a count vector of k-mers in the two strings
#'
#' @param str1 a character vector specifying the first string to be tested
#' @param str2 a character vector specifying the second string to be tested
#' @param k size of k-mer
#'
#' @return a list containing two count vectors
#'
#' @examples
#' alfa::countVectors("ATGTGAC", "ATGACTT", 3)
#'
#' @export
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
#' @return the Euclidean distance between str1 and str2
#'
#' @examples
#' alfa::euclideanDistance("ATCATC", "ATTATC", 3)
#'
#' @references
#' Zielezinski, A., Vinga, S., Almeida, J. et al. Alignment-free sequence
#' comparison: benefits, applications, and tools. Genome Biol 18, 186 (2017).
#'
#' @export
euclideanDistance <- function(str1, str2, k) {
  vecs <- countVectors(str1, str2, k)
  stopifnot(length(vecs) == 2)
  str1Count <- vecs$first
  str2Count <- vecs$second
  distMat <- dist(rbind(str1Count, str2Count))
  return(distMat[1])
}

#' Calculate the Euclidean distance between two strings
#' using pre-computed values
#'
#' @param countVec1 the precomputed count vector for str1
#' @param countVec2 the precomputed count vector for str2
#'
#' @return the Euclidean distance between countVec1 and countVec2
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
#' @return the standardized Euclidean distance between str1 and str2
#'
#' @examples
#' alfa::standardizedEuclidean("ACAATC", "ACATAA", 3)
#'
#' @references
#' Wu TJ, Burke JP, Davison DB. A measure of DNA sequence dissimilarity
#' based on Mahalanobis distance between frequencies of words.
#' Biometrics. 1997 Dec;53(4):1431-9. PMID: 9423258.
#'
#' @export
standardizedEuclidean <- function(str1, str2, k) {
  vecs <- standardizedFreq(str1, str2, k)
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
#' @return the standardized Euclidean distance between str1 and str2
fastStdEuclideanDistance <- function(countVec1, countVec2, str1Len, str2Len,
                                     words, overlapCap, k) {
  freqs1 <- countVec1 / sum(countVec1)
  freqs2 <- countVec2 / sum(countVec2)
  # var1 <- vector(length = length(freqs1))
  # var2 <- vector(length = length(freqs2))
  var1 <- lapply(overlapCap, freqVariance, seqLen = str1Len, wordLen = k)
  var2 <- lapply(overlapCap, freqVariance, seqLen = str2Len, wordLen = k)
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
#' @return the Hamming distance between str1 and str2
hammingDistance <- function(str1, str2, k) {
  vecs <- countVectors(str1, str2, k)
  stopifnot(length(vecs) == 2)
  str1Count <- vecs$first
  str2Count <- vecs$second
  return(mean(str1Count != str2Count))
}
