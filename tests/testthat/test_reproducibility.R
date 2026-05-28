test_that("results are reproducible", {
  set.seed(123)
  X <- list(
    matrix(rnorm(250), 10, 25),
    matrix(rnorm(120), 10, 12)
  )
  Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
  end <- c(0.3, 0.3, 3)
  step <- c(0.1, 0.1, 0.5)
  
  set.seed(3)
  res1 <- getCCA3(X, Z, end, step, numCV = 4)
  
  set.seed(3)
  res2 <- getCCA3(X, Z, end, step, numCV = 4)
  
  expect_equal(res1$corr, res2$corr)
  for (i in seq_along(res1$cc3.weight.x)) {
    
    expect_equal(
      as.vector(res1$cc3.weight.x[[i]]),
          as.vector(res2$cc3.weight.x[[i]]))
    
    expect_equal(
      as.vector(res1$cc3.CV.x[[i]]),
          as.vector(res2$cc3.CV.x[[i]]))
    
  }
  expect_equal(
    as.vector(res1$cc3.weight.z),
        as.vector(res2$cc3.weight.z))  
    expect_equal(
      as.vector(res1$cc3.CV.z),
          as.vector(res2$cc3.CV.z)) 
    
  
})