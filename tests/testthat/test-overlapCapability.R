test_that("overlapCapability works", {
  expect_equal(alfa::overlapCapability("ACAC"), c(0,1,0,1))
  expect_equal(alfa::overlapCapability("AAAC"), c(0,0,0,1))
  expect_equal(alfa::overlapCapability("ATGC"), c(0,0,0,1))
})
