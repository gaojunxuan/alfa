test_that("allKMers works", {
  kmers <- alfa::allKMers("ATGTGAC", 3)
  expect_equal(kmers,
               c("ATG", "TGT", "GTG", "TGA", "GAC"))
})

test_that("allKMers fails", {
  expect_error(alfa::allKMers("ATGTGAC", 30))
})
