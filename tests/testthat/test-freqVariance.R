test_that("freqVariance works", {
  expect_equal(round(alfa:::freqVariance(3, 3, c(1,1,1)), 3),
               round(0.015869140625, 3))
})
