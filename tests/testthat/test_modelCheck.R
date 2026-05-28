#foreach::registerDoSEQ()

## Output structure test for two biological datasets
test_that("works with minimum number of blocks", {
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12)
  )
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  res <- getCCA3(X, Z, end, step, numCV = 4)
  
  
  expect_true(is.list(res))
  expect_true(length(res$corr) == 4)
})

## Output structure test for three biological datasets

test_that("works with 3 datasets", {
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12),
    matrix(rnorm(150), 10, 15)
  )
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.1, 0.5)
  
  res <- getCCA3(X, Z, end, step, numCV = 4)
  
  expect_true(is.list(res))
  expect_true(length(res$corr) == 4)
})