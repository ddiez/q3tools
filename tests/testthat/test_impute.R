context("Impute return value")
test_that("imputeGroup returns matrix", {
  expect_is(imputeGroup(matrix(NA, ncol = 2, nrow = 2), c(1, 2)), "matrix")
  expect_false(any(is.na(imputeGroup(matrix(NA, ncol = 2, nrow = 2), c(1, 2)))))
})
