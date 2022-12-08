test_that("countVectors works", {
  cntVec <- alfa::countVectors("ATGTGAC", "ATGACTT", 3)
  expect_equal(sum(alfa::countVectors("ATGTGAC", "ATGACTT", 3)$first), 5)
  expect_equal(sum(alfa::countVectors("ATGTGAC", "ATGACTT", 3)$second), 5)
  expect_equal(alfa::countVectors("ATGTGAC", "ATGACTT", 3)$first[[3]], 1)
})
