#' Calculate supervised CCA on three data matrices
#'
#' getCCA3() starts the supervised sparse CCA 
#' 
#' 
#' 
#' 
#' The getCCA3() function performs sparse 
#' canonical correlation analysis on three data matrices 
#' (two biological data sets X,Y, and one design data set Z) with elastic net. 
#' Data sets X and Y are strongly regularized (ridge regression); Z if neccessary (det(Z)==0).
#' 
#' It normalizes the data matrices, computes canonical variables, 
#' and updates the data matrices by removing the latent variables iteratively 
#' for a specified number of canonical variables (numCV). 
#' It returns a list containing matrices with weight vectors and canonical variables, 
#' correlation coefficients, and optimal sparsity parameters (lambda values) 
#' for each canonical variable.
#' 
#' @param X biological data set with dim(X): n rows, p features
#' @param Y biological data set with dim(Y): n rows, q features
#' @param Z design data set with dim(Z): n rows, r features; describes experimental design with binary vectors
#' @param end,step end and stepsize for lambda.{x,y,z} grid search.
#' @param n.r number of resampling runs, 10 is okay
#' @param max.counter.test number of random start vectors for iteration (inner loop)
#' @param numCV number of canonical variables which will be calculated
#' @param alpha control the regularization strength 
#'
#' @return A list with the elements 
#'   cc3.weights.{xyz}: Matrices with weights of canonical variables for data set X, Y and Z respectively (columnwise)
#'   cc3.CV.{xyz}: matrices with canonical variables for data set X, Y and Z respectively (columnwise)
#'   corr: matrix with correlation coefficient for pairwise (X-Y,X-Z,YZ) and correlation of can. variable: sum(corr(X-Z),corr(Y-Z))/2 (columnwise)
#'         STN: contradictory documentation, corr could be: Vector with absolute correlation coefficients of canonical variables (columnwise): (corr(cv(X),cv(Z))+ corr(cv(Y),cv(Z))/2    
#'   lambda: best lambda.{xyz} for each data set and for each can. var. (columnwise)
#'   numCV: Number of canonical variables
#'   
#' @keywords supervised CCA sparse
#' @author Andrea Thum
#' @importFrom MASS ginv
#' @export
#'
#' @examples
#' TRUE
getCCA3 <- function(X, Y, Z,
                    end = c(0.3, 0.5, 3),
                    step = c(0.01, 0.01, 0.2),
                    numCV = 10,
                    n.r = 10,
                    max.counter.test = 10,
                    alpha=1) {
  if (dim(X)[1] != dim(Y)[1] || dim(Y)[1] != dim(Z)[1])
    stop("X,Y,Z data matrices have different sample sizes: ", 
         dim(X)[1], " ", dim(Y)[1], " ", dim(Z)[1])
  
  cc3.weight.x <- matrix(ncol = 0, nrow = dim(X)[2])
  cc3.weight.y <- matrix(ncol = 0, nrow = dim(Y)[2])
  cc3.weight.z <- matrix(ncol = 0, nrow = dim(Z)[2])
  cc3.CV.x <- matrix(ncol = 0, nrow = dim(X)[1])
  cc3.CV.y <- matrix(ncol = 0, nrow = dim(X)[1])
  cc3.CV.z <- matrix(ncol = 0, nrow = dim(X)[1])
  
  rownames(cc3.weight.x) <- colnames(X)
  rownames(cc3.weight.y) <- colnames(Y)
  rownames(cc3.weight.z) <- colnames(Z)
  rownames(cc3.CV.x) <- rownames(X)
  rownames(cc3.CV.y) <- rownames(X)
  rownames(cc3.CV.z) <- rownames(X)
  
  all.lambdas <- matrix(ncol = 0, nrow = 3)
  all.corr <- c()
  
  Z <- apply(Z, 2, function(y) {
    if (var(y) == 0)
      y - mean(y)
    else
      (y - mean(y)) / var(y)
  })
  
  dims <- c(dim(X)[2], dim(Y)[2], dim(Z)[2])
  
  # determine pseudo matrix once - Z does not change
  if (abs(det(var(Z))) > 10 ^ -20) {
    # if var(Z) is invertible
    #Zp <- solve(var(Z)) %*% t(Z)
    Zp <- ginv(var(Z)) %*% t(Z)
    
    } else{
    Zp <- (diag(1 / sqrt(diag(var(Z)))) * alpha) %*% t(Z) # otherwise: regularize
  }
  
  canVar = 1
  while (canVar <= numCV) {
    # calculate the canonical variables
    # normalize data matrices
    X <- apply(X, 2, function(y) {
      if (var(y) == 0)
        y - mean(y)
      else
        (y - mean(y)) / var(y)})
    
    Y <- apply(Y, 2, function(y) {
      if (var(y) == 0)
        y - mean(y)
      else
        (y - mean(y)) / var(y)})
    
    # determine Pseudomatrices for every canonical correlation step, 
    # including ridge regression
    Xp <- (diag(1 / sqrt(diag(var(X)))) * alpha) %*% t(X) # strong regularization
    Yp <- (diag(1 / sqrt(diag(var(Y)))) * alpha) %*% t(Y)
    
    # Pseudomatrices * ...
    XpZ <- Xp %*% Z
    YpZ <- Yp %*% Z
    ZpX <- Zp %*% X
    ZpY <- Zp %*% Y
    
    # get best combination of sparsity parameters
    results <- get.best.lambdas(X, Y, Z,
      end = end,
      n.r = n.r,
      step = step,
      max.counter.test = max.counter.test) 
    
    if (is.null(results)) {
      break
    }
    lambdas <- c(results$best.lambda.x,
                 results$best.lambda.y,
                 results$best.lambda.z)
    corr.scca <- results$corr
    
    xj <- results$bestVector[1:dims[1]]
    yj <- results$bestVector[(dims[1] + 1):(dims[1] + dims[2])]
    zj <- results$bestVector[(dims[1] + dims[2] + 1):(dims[1] + dims[2] + dims[3])]
    
    # update data matrices by removing the latent variable - only for X and Y
    zi <- Z %*% zj
    xi <- X %*% xj
    reg <- apply(X, 2, function(x) {lm(x ~ xi)} )
    X <- sapply(reg, function(x) {x[[2]]} )
    
    yi <- Y %*% yj
    reg <- apply(Y, 2, function(x) {lm(x ~ yi)} )
    Y <- sapply(reg, function(x) {x[[2]]} )
    
    all.lambdas <- cbind(all.lambdas, lambdas)
    all.corr <- c(all.corr, corr.scca)
    cc3.weight.x <- cbind(cc3.weight.x, xj)
    cc3.weight.y <- cbind(cc3.weight.y, yj)
    cc3.weight.z <- cbind(cc3.weight.z, zj)
    cc3.CV.x <- cbind(cc3.CV.x, xi)
    cc3.CV.y <- cbind(cc3.CV.y, yi)
    cc3.CV.z <- cbind(cc3.CV.z, zi)
    
    canVar <- canVar + 1
  } # while canVar
  
  return(list(cc3.weight.x = cc3.weight.x,
              cc3.weight.y = cc3.weight.y,
              cc3.weight.z = cc3.weight.z,
              cc3.CV.x = cc3.CV.x,
              cc3.CV.y = cc3.CV.y,
              cc3.CV.z = cc3.CV.z,
              corr = all.corr,
              lambda = all.lambdas,
              num.CV = numCV))
}
