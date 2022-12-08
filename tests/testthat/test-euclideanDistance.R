test_that("euclideanDistance works", {
  expect_equal(round(alfa::euclideanDistance("ATCATC", "ATTATC", 3), 5), 2.44949)
})

test_that("euclideanDistance fails", {
  expect_error(alfa::euclideanDistance("ATCATC", "ATTATC", 10))
})
