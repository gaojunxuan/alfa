test_that("readFASTA fails on non-existent file", {
  expect_error(alfa::readFASTA("aasdjfalsdk"))
})
