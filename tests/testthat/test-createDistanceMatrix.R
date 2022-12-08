test_that("createDistanceMatrix works", {
  testStringSet <- Biostrings::DNAStringSet(c("AAA","ATT","ATC"))
  distMat <- alfa::createDistanceMatrix(testStringSet, k=3)
  expected <- matrix(NA, ncol = 3, nrow = 3)
  expected[lower.tri(expected)] <- 1.414214
  expected[upper.tri(expected)] <- 1.414214
  expect_equal(round(distMat, 6),
               expected)
})

test_that("createDistanceMatrix fails", {
  testStringSet <- Biostrings::DNAStringSet(c("AAA","ATT","ATC"))
  expect_error(alfa::createDistanceMatrix(testStringSet, k=10))
})
