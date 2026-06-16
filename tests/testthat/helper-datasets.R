set.seed(2)
X <- list(
  matrix(rnorm(250), 10, 25),
  matrix(rnorm(120), 10, 12)
)
Z <- matrix(sample(c(0, 1), 50, replace = TRUE), ncol = 5)  
end <- c(0.3, 0.3, 3)
step <- c(0.1, 0.1, 0.5)