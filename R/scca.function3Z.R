#' Title scca.function3Z
#'
#' @param XpZ,YpZ,ZpX,ZpY XpZ = regularized Pseudoinverse(X) * Z, analogously: YpZ, ZpX, ZpY with X,Y: biol. data sets; Z: design data set
#' @param dims 
#' @param lambda.x,lambda.y,lambda.z sparsity parameters for Data Sets X, Y, Z
#' @param max.iter max. number of iterations for power method
#' @param shift lasso-regularization after {shift} iterations
#' @param x.initial,y.initial,z.initial may start with preset start vectors, if NULL: random
#'
#' @return a list with the components:
#'   {xyz}.new: eigenvectors, i.e. weight vectors for supervised sparse CCA for 3 data sets (X,Y;Z) and given sparsity (lambdas)
#    i: number of iterations (should not be max.iter!)
#    null: at least one sparsity parameter lambda is too strict: weight vector = (0,...,0)

#' @importFrom stats cor lm median runif var
#' @export
#'
#' @examples
#' TRUE

scca.function3Z <- function(XpZ, YpZ, ZpX, ZpY,
                            dims,
                            lambda.x, lambda.y, lambda.z,
                            max.iter = 500,
                            shift = 2,
                            x.initial = NULL, y.initial = NULL, z.initial = NULL)
{
  p <- dims[1]
  q <- dims[2]
  r <- dims[3]
  
  eps <- 0.0001	# convergence criterion
  
  isnull = FALSE
  
  # random start vectors if init vectors are unknown
  if (is.null(z.initial))
    z.initial <- runif(r)
  if (is.null(x.initial))
    x.initial <- runif(p)
  if (is.null(y.initial))
    y.initial <- runif(q)
  
  #normalize
  x.initial <- x.initial / as.numeric(sqrt(t(x.initial) %*% x.initial))
  y.initial <- y.initial / as.numeric(sqrt(t(y.initial) %*% y.initial))
  z.initial <- z.initial / as.numeric(sqrt(t(z.initial) %*% z.initial))
  
  diff.x <- eps * 10
  diff.y <- eps * 10
  diff.z <- eps * 10
  
  nextstart = FALSE
  i <- 0		# number of iterations used by SCCA
  
  while ((i < max.iter) &
         ((diff.x > eps) || (diff.y > eps) || (diff.z > eps)))
  {
    i <- i + 1
    # Update x-vector
    
    hx <- x.initial + XpZ %*% z.initial #update
    length.hx <- as.numeric(sqrt(t(hx) %*% hx))
    if (is.na(length.hx) |
        length.hx == 0) {
      isnull = TRUE
      break
    }
    hx <- hx / length.hx
    x.new <- hx
    if (i > shift) {
      x.new <- abs(hx) - 0.5 * lambda.x #lasso
      x.new <- (x.new + abs(x.new)) / 2
      x.new <- x.new * sign(hx)
      length.x.new <- as.numeric(sqrt(t(x.new) %*% x.new))
      if (is.na(length.hx) |
          length.x.new == 0) {
        isnull = TRUE
        break
      }
      x.new <- x.new / length.x.new
    }
    
    # Update y-vector
    
    hy <- y.initial + YpZ %*% z.initial
    length.hy <- as.numeric(sqrt(t(hy) %*% hy))
    if (is.na(length.hy) |
        length.hy == 0) {
      isnull = TRUE
      break
    } # print(sqrt(t(hy)%*%hy));}
    hy <- hy / length.hy
    y.new <- hy
    if (i > shift) {
      y.new <- abs(hy) - 0.5 * lambda.y
      y.new <- (y.new + abs(y.new)) / 2
      y.new <- y.new * sign(hy)
      length.y.new <- as.numeric(sqrt(t(y.new) %*% y.new))
      if (is.na(length.hy) |
          length.y.new == 0) {
        isnull = TRUE
        break
      }
      y.new <- y.new / length.y.new
    }
    
    # Update z-vector, Z: design data set
    
    # here we need the new values x.new and y.new to prevent switching 
    # between two eigen vectors
    hz <- ZpX %*% x.new + ZpY %*% y.new +  z.initial
    length.hz <- as.numeric(sqrt(t(hz) %*% hz))
    if (is.na(length.hz) |
        length.hz == 0) {
      isnull = TRUE
      break
    }
    hz <- hz / length.hz
    z.new <- hz
    if (i > shift) {
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
    diff.x <- max(abs(x.initial - x.new))
    diff.y <- max(abs(y.initial - y.new))
    diff.z <- max(abs(z.initial - z.new))
    
    x.initial <- x.new
    y.initial <- y.new
    z.initial <- z.new
    
    if ((max(abs(z.new)) == 0) ||
        (max(abs(x.new)) == 0) ||
        (max(abs(y.new)) == 0)) {
      isnull = TRUE
      break
    }
    
  } # while diff.x
  
  list(x.new = x.initial,
      y.new = y.initial,
      z.new = z.initial,
      i = i,
      null = isnull)
}
