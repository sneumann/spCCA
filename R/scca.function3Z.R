#' Title scca.function3Z
#' 
#' Implements the modified iterative power method for three (regularized) data sets 
#' with lasso penalty for given lasso parameters lambda.(x,y,z).
#' Pseudomatrices are used as input for faster calculation.
#' This function is called by get.best.lambdas()
#'
#' @param XpZ,ZpX For regularized pseudomatrix Xp of X: regularize(t(X)*X)) * t(X); then XpZ = Xp * Z is input. Z is the design data set.
#' @param dims number of features of all data sets including Z
#' @param lambda.x,lambda.z lasso parameters lambda for data sets. lambda = 0: no sparsity
#' @param max.iter max. number of iterations for power method
#' @param shift start lasso regularization after shift iteration steps
#' @param x.initial,z.initial initial vectors for iteration. If NULL: filled with random values
#'
#' @return a list with the components:
#'   (xz).new: eigenvectors, i.e. weight vectors for supervised sparse CCA for all data sets (X;Z) and given sparsity (lambdas)
#    n: number of iterations until convergence (should not be max.iter!)
#    null: true if at least one lambda is too high and all weights of one weight vector are zero}, i.e. at least one sparsity parameter lambda is too strict: weight vector = (0,...,0)
#'
#' @note Test output variable n, if iteration did converge (n<max.iter)!
#' @importFrom stats cor lm median runif var
#' @importFrom utils tail
#' @export
#' @author Andrea Thum, Elena Parkhomenko
#' @keywords CCA lasso
#'
#' @examples
#' TRUE

scca.function3Z <- function(XpZ, ZpX,
                            dims,
                            lambda.x, lambda.z,
                            max.iter = 100,
                            shift = 2,
                            x.initial = list(), z.initial = NULL)
{
  
  eps <- 0.0001	# convergence criterion
  
  isnull = FALSE
  
  # random start vectors if init vectors are unknown
  if (is.null(z.initial))
    z.initial <- runif(tail(dims,1))
  z.initial <- z.initial / as.numeric(sqrt(t(z.initial) %*% z.initial))  # normalize
  diff.z <- eps * 10
  
  if (is.list(x.initial) && length(x.initial) == 0){
    for (i in seq_len(length(dims)-1)){
      x.initial[[i]] <- runif(dims[i])
    }
  }
  
  diff.x <- list()
  for (i in seq_along(x.initial)){
    x.initial[[i]] <- x.initial[[i]] / as.numeric(sqrt(t(x.initial[[i]]) %*% x.initial[[i]])) # normalize
    diff.x[i] <- eps * 10
  }
  
  nextstart = FALSE
  n <- 0	# number of iterations used by SCCA
  x.new <- list()
  z.new <- list()
  
  while ((n < max.iter) & (any(c(diff.x, diff.z) > eps)))
  {
    n <- n + 1
    # Update x-vector
    for (i in seq_along(x.initial)){
      hx <- x.initial[[i]] + XpZ[[i]] %*% z.initial #update
      length.hx <- as.numeric(sqrt(t(hx) %*% hx))
      if (is.na(length.hx) |
          length.hx == 0) {
        isnull = TRUE
        break
      }
      hx <- hx / length.hx
      x.new[[i]] <- hx
      if (n > shift) {
        x.new[[i]] <- abs(hx) - 0.5 * lambda.x[i] #lasso
        x.new[[i]] <- (x.new[[i]] + abs(x.new[[i]])) / 2
        x.new[[i]] <- x.new[[i]] * sign(hx)
        length.x.new <- as.numeric(sqrt(t(x.new[[i]]) %*% x.new[[i]]))
        if (is.na(length.hx) |
            length.x.new == 0) {
          isnull = TRUE
          break
        }
        x.new[[i]] <- x.new[[i]] / length.x.new
      }
    }
    
    # Update z-vector, Z: design data set
    
    # here we need the new values x.new and y.new to prevent switching 
    # between two eigen vectors
    hz <- matrix(0, nrow = nrow(ZpX[[1]]), ncol = 1)
    for (i in seq_along(x.new)){
      hz <- hz + ZpX[[i]] %*% x.new[[i]]
    }
    hz <- hz +  z.initial
    length.hz <- as.numeric(sqrt(t(hz) %*% hz))
    if (is.na(length.hz) | length.hz == 0) {
      isnull = TRUE
      break
    }
    hz <- hz / length.hz
    z.new <- hz
    if (n > shift) {
      z.new <- abs(hz) - 0.5 * lambda.z
      z.new <- (z.new + abs(z.new)) / 2
      z.new <- z.new * sign(hz)
      length.z.new <- as.numeric(sqrt(t(z.new) %*% z.new))
      if (is.na(length.hz) |
          length.z.new == 0) {
        isnull = TRUE
        break
      }
      z.new <- z.new / length.z.new
    }
    
    # Convergence measures
    for (i in seq_along(x.initial)){
      diff.x[[i]] <- max(abs(x.initial[[i]] - x.new[[i]]))
    }
    diff.z <- max(abs(z.initial - z.new))
    
    x.initial <- x.new
    z.initial <- z.new
    
    if ((max(abs(z.new)) == 0) ||  (any(sapply(x.new, function(x) (max(abs(x)) == 0)))))  { 
      isnull = TRUE
      break
    }
    
  } # while diff.x
  
  return(list(x.new = x.initial,
              z.new = z.initial,
              n = n,
              null = isnull))
}
