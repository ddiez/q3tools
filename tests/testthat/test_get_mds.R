context("Obtain MDS values")
test_that("get_mds returns MDS object", {
  expect_is(get_mds(mtcars), "MDS")
})
