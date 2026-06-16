
## 1. Output structure test for two biological datasets
test_that("works with minimum number of blocks", {
  
  res <- getCCA3(X, Z, end, step, numCV = 4)
  
  expect_true(is.list(res))
  expect_true(length(res$corr) == 4)
})

## 2. Output structure test for three biological datasets

test_that("works with 3 datasets", {
  X3 <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12),
    matrix(rnorm(150), 10, 15)
  )
  set.seed(123)
  Z3 <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.1, 0.5)
  
  res <- getCCA3(X3, Z3, end, step, numCV = 4)
  
  expect_true(is.list(res))
  expect_true(length(res$corr) == 4)
})

## 3. More features than samples: common in omics data
test_that("high-dimensional setting runs", {
  
  X <- list(
    matrix(rnorm(1000), 10, 100),
    matrix(rnorm(500), 10, 50)
  )
  
  Z <- matrix(sample(0:1, 50, TRUE), 10, 5)
  
  expect_no_error(
    getCCA3(X, Z, end, step, numCV = 2)
  )
})

## 4. Swap-order invariance test
test_that("results are invariant to ordering of X blocks", {
  
  X1 <- matrix(rnorm(250), 10, 25)
  X2 <- matrix(rnorm(120), 10, 12)
  
  set.seed(1)
  res1 <- getCCA3(list(X1,X2), Z, end, step, numCV = 4)
  set.seed(1)
  res2 <- getCCA3(list(X2,X1), Z, end, step, numCV = 4)
  
  # compare feature selection proportions for the first CV
  
  wA <- res1$cc3.weight.x[[1]][,1]
  wB <- res2$cc3.weight.x[[2]][,1]  # swapped
  selA <- which(wA != 0)
  selB <- which(wB != 0)
  
  # Jaccard overlap
  jaccard <- length(intersect(selA, selB)) / length(union(selA, selB))
  
  # Weight similarity
  weight_cor <- abs(cor(wA, wB))
  
  # Sparsity proportion
  sparsity_diff <- abs(mean(wA != 0) - mean(wB != 0))
  
  expect_true(jaccard > 0.8)
  expect_true(weight_cor > 0.8)
  expect_true(sparsity_diff < 0.1)
})

## 5. Edge cases
test_that("handles single feature case", {
  
  X_1 <- list(
    matrix(rnorm(10), 10, 1),
    matrix(rnorm(10), 10, 1)
  )
  set.seed(123)
  Z_1 <- matrix(sample(c(0, 1), 10, replace = TRUE), ncol = 1)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  res <- getCCA3(X_1, Z_1, end, step, numCV = 1)
  expect_length(res$corr, 1)
})

## 6. Extremely sparse lambda
test_that("large penalties do not crash", {
  
  end <- c(10,10,10)
  step <- c(10,10,10)
  
  expect_no_error(
    getCCA3(X, Z, end, step, numCV=2)
  )
})

## 7. Tiny penalties
test_that("small penalties do not crash", {
  
  end <- c(0.001,0.001,0.001)
  step <- c(0.001,0.001,0.001)
  
  expect_no_error(
    getCCA3(X, Z, end, step, numCV=2)
  )
})

## 8. Constant feature: variance = 0 is common source of NaN correlations
test_that("constant feature handled gracefully", {
  
  X <- list(
    cbind(rep(1, 10), matrix(rnorm(40), 10, 4)),
    matrix(rnorm(50), 10, 5)
  )
  
  Z <- matrix(sample(0:1, 50, TRUE), 10, 5)
  
  expect_no_error(
    getCCA3(X, Z, end, step, numCV=2)
  )
})