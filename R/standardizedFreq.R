overlapCapability <- function(word) {
  res <- c()
  wordLen <- stringr::str_length(word)
  for (i in 1:wordLen) {
    v <- 0
    if (substr(word, 1, i) == substr(word, wordLen - i + 1, wordLen))
      v <- 1
    res <- c(res, v)
  }
  return(res)
}

freqVariance <- function(seqLen, wordLen, overlap) {
  wordProbability <- c()
  currProb <- 1
  for (i in 1:wordLen) {
    wordProbability <- c(wordProbability, 0.25 * currProb)
    currProb <- currProb * 0.25
  }
  p <- wordProbability[length(wordProbability)]
  sumRange <- 1:min(seqLen - wordLen, wordLen - 1)
  summation <- c()
  for (k in sumRange) {
    val <- (seqLen - wordLen + 1 - k) * overlap[wordLen - k] *
      wordProbability[k]
    summation <- c(summation, val)
  }
  sumVal <- sum(summation)
  firstTerm <- p * (seqLen - wordLen + 1) * (1 - p * (seqLen - wordLen + 1))
  secondTerm <- p ^ 2 * (seqLen - wordLen + 1 - wordLen) *
    (seqLen - wordLen + 2 - wordLen)
  thirdTerm <- 2 * p * sumVal
  return(firstTerm + secondTerm + thirdTerm)
}

standardizedFreq <- function(str1, str2, k) {
  str1Len <- stringr::str_length(str1)
  str2Len <- stringr::str_length(str2)
  kMerCounts <- alfa::countVectors(str1, str2, k)
  freqs1 <- kMerCounts$first / sum(kMerCounts$first)
  freqs2 <- kMerCounts$second / sum(kMerCounts$first)
  overlaps1 <- c()
  overlaps2 <- c()
  wordTableDf1 <- as.data.frame(freqs1)
  wordTableDf2 <- as.data.frame(freqs2)
  for (word in wordTableDf1$Var1) {
    overlaps1 <- cbind(overlaps1, overlapCapability(word))
  }
  for (word in wordTableDf2$Var1) {
    overlaps2 <- cbind(overlaps2, overlapCapability(word))
  }
  var1 <- c()
  var2 <- c()
  for (i in 1:length(freqs1)) {
    var <- freqVariance(str1Len, stringr::str_length(wordTableDf1$Var1[i]), overlaps1[,i])
    var1 <- c(var1, var)
  }
  for (i in 1:length(freqs2)) {
    var <- freqVariance(str2Len, stringr::str_length(wordTableDf2$Var1[i]), overlaps2[,i])
    var2 <- c(var2, var)
  }
  return(list("first" = freqs1 / var1, "second" = freqs2 / var2))
}
