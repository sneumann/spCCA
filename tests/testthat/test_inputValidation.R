## Input validation test

test_that("Errors when datasets in X have different row counts", {
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 12, 10),
    matrix(rnorm(150), 10, 15)
  )
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.01, 0.01, 0.2)
  
  expect_error(getCCA3(X, Z, end, step, numCV = 4))
  
  
})

test_that("Errors when X and Z row counts differ", {
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12),
    matrix(rnorm(150), 10, 15)
  )
  Z <- matrix(sample(c(0, 1), 10, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.01, 0.01, 0.2)
  
  expect_error(getCCA3(X, Z, end, step, numCV = 4))
  
  
})