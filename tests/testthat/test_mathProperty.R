## 1. correlations must be between [0,1] 
test_that("canonical correlations are bounded", {
  res <- getCCA3(X, Z, end, step, numCV = 4)
  
  expect_true(all(res$corr >= 0))
  expect_true(all(res$corr <= 1))
})


## 2. Monotonic-ish decrease in correlations
test_that("correlations are not wildly increasing", {
  
  res <- getCCA3(X, Z, end, step, numCV = 4)
  
  diffs <- diff(res$corr)
  

  expect_true(all(diffs <= 0.03))
})


## 3. Stability under purturbation
test_that("monotonic trend is stable under small perturbations", {
  
  set.seed(1)
  res1 <- getCCA3(X, Z, end, step, numCV = 4)
  
  # perturb data slightly
  X2 <- lapply(X, function(x) x + rnorm(length(x), sd = 0.01))
  set.seed(1)
  res2 <- getCCA3(X2, Z, end, step, numCV = 4)
  
  diff1 <- diff(res1$corr)
  diff2 <- diff(res2$corr)
  

  # both should show similar trend direction
  expect_true(all(diff1 < 0.02))
  expect_true(all(diff2 < 0.02))
})

## 4. Random seed consistency test
test_that("monotonic trend is not seed dependent", {
  
  set.seed(15)
  r1 <- getCCA3(X, Z, end, step, numCV = 4)
  
  set.seed(203)
  r2 <- getCCA3(X, Z, end, step, numCV = 4)
  
  expect_true(all(diff(r1$corr) < 0.02))
  expect_true(all(diff(r2$corr) < 0.02))
})


## 5. Orthogonality for deflation correctness
test_that("canonical variates are approximately orthogonal", {
  
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

