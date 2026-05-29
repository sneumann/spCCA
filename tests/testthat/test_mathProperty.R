## 1. correlations must be between [0,1] 
test_that("canonical correlations are bounded", {
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12)
  )
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  res <- getCCA3(X, Z, end, step, numCV = 4)
  
  expect_true(all(res$correlations >= 0))
  expect_true(all(res$correlations <= 1))
})


## 2. Monotonic-ish decrease in correlations
test_that("correlations are not wildly increasing", {
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12)
  )
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  res <- getCCA3(X, Z, end, step, numCV = 4)
  
  diffs <- diff(res$corr)

  expect_true(all(diffs <= 0.015))
})


## 3. Stability under purturbation
test_that("monotonic trend is stable under small perturbations", {
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12)
  )
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  set.seed(1)
  res1 <- getCCA3(X, Z, end, step, numCV = 4)
  
  # perturb data slightly
  X2 <- lapply(X, function(x) x + rnorm(length(x), sd = 0.01))
  set.seed(1)
  res2 <- getCCA3(X2, Z, end, step, numCV = 4)
  
  diff1 <- diff(res1$corr)
  diff2 <- diff(res2$corr)
  

  # both should show similar trend direction
  expect_true(all(diff1 < 0.015))
  expect_true(all(diff2 < 0.015))
})

## 4. Random seed consistency test
test_that("monotonic trend is not seed dependent", {
  
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12)
  )
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  set.seed(15)
  r1 <- getCCA3(X, Z, end, step, numCV = 4)
  
  set.seed(203)
  r2 <- getCCA3(X, Z, end, step, numCV = 4)
  
  expect_true(all(diff(r1$corr) < 0.02))
  expect_true(all(diff(r2$corr) < 0.02))
  expect_true(cor(r1$corr, r2$corr) > 0.7)
})


## 5. Orthogonality for deflation correctness
test_that("canonical variates are approximately orthogonal", {
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12)
  )
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  res <- getCCA3(X, Z, end, step, numCV = 4)
  
  for (i in seq_along(res$cc3.CV.x)) {
    
    CV <- res$cc3.CV.x[[i]]  # n samples ?? numCV matrix
    
    corr_mat <- cor(CV)
    
    off_diag <- corr_mat[lower.tri(corr_mat)]
    
    expect_true(all(abs(off_diag) < 0.1e-10))  # must show low cross component correlation
  }
  
  CVz <- res$cc3.CV.z
  corr_z <- cor(CVz)
  off_diag_z <- corr_z[lower.tri(corr_z)]

  expect_true(mean(abs(off_diag_z)) < 0.5)  ## Expecting strong orthogonality in Z is unrealistic, since design matrix naturally have lower-rank structure.
  
})



## 6. scaling invariance test
test_that("results are stable under scaling of X list", {
  
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12)
  )
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  set.seed(1)
  res1 <- getCCA3(X, Z, end, step, numCV = 4)
  
  X_scaled <- lapply(X, function(x) x * 10)
  set.seed(1)
  res2 <- getCCA3(X_scaled, Z, end, step, numCV = 4)
  
  X_scaled2 <- lapply(X, function(x) x * 2)
  set.seed(1)
  res3 <- getCCA3(X_scaled2, Z, end, step, numCV = 4)

  expect_true(
    all.equal(res1, res2, tolerance = 1e-7)
  )
  expect_true(
    all.equal(res1, res3, tolerance = 1e-7)
  )
})

## 7. Swap-order invariance test
test_that("results are invariant to ordering of X blocks", {
  
  X1 <- matrix(rnorm(250), 10, 25)
  X2 <- matrix(rnorm(120), 10, 12)
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
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

# 8. Edge cases
test_that("handles single feature case", {
  
  X <- list(
    matrix(rnorm(10), 10, 1),
    matrix(rnorm(10), 10, 1)
  )
  set.seed(123)
  Z <- matrix(sample(c(0, 1), 10, replace = TRUE), ncol = 1)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  res <- getCCA3(X, Z, end, step, numCV = 1)
  expect_length(res$corr, 1)
})
  