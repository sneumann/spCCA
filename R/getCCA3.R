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
#' The getCCA3 first normalizes the data matrices (column-wise mean centering and unit variance), 
#' computes canonical variables, 
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
getCCA3 <- function(X, Z,
                    end, # c(0.3, 0.5, 3),
                    step, # c(0.01, 0.01, 0.2),
                    numCV = 10,
                    n.r = 10,
                    max.counter.test = 10) {
  if (length(end) != length(X)+1 || length(step) != length(X)+1)
    stop('Length of end or stepsize does not match the total number of datasets')
  
  Total_data <- c(X,list(Z))
  if (length(unique(sapply(Total_data,nrow)))!=1)
    stop('Data Matrices have different sample sizes: ',sapply(Total_data,nrow))
  n_sets <- length(X)
  
  cc3.weight <- vector("list", n_sets)
  cc3.CV <- vector("list", n_sets)
  
  for (i in seq_len(n_sets)) { 
    cc3.weight[[i]] <- matrix(ncol = 0, nrow = ncol(X[[i]]))
    rownames(cc3.weight[[i]]) <- colnames(X[[i]])
    
    cc3.CV[[i]] <- matrix(ncol = 0, nrow = nrow(X[[i]]))
    rownames(cc3.CV[[i]]) <- rownames(X[[i]])
  }
  
  cc3.weight.z <- matrix(ncol = 0, nrow = dim(Z)[2])
  cc3.CV.z <- matrix(ncol = 0, nrow = dim(Z)[1])
  rownames(cc3.weight.z) <- colnames(Z)
  rownames(cc3.CV.z) <- rownames(Z)
  
  all.lambdas <- matrix(ncol = 0, nrow = length(Total_data))
  all.corr <- c()
  
  Z_std <- apply(Z, 2, function(y) {
    if (sd(y) == 0)
      y - mean(y)
    else
      (y - mean(y)) / sd(y)
  })
  
  dims <- sapply(Total_data, ncol)
  
  canVar = 1
  while (canVar <= numCV) {
   
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
                 results$best.lambda.z)
    corr.scca <- results$corr
    
    # split vectors
    split_idx <- rep(seq_along(dims), dims)
    vectors <- split(results$bestVector, split_idx)
    
    zj <- vectors[[length(vectors)]]
    
    xj_list <- vectors[-length(vectors)]
    
    # standardize X before deflation
    standardize_X <- function(M){
      mu <- colMeans(M)
      sdv <- apply(M, 2, sd)
      sdv[sdv == 0] <- 1
      
      M_scaled <- sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
      M_scaled
    }
    for (i in seq_along(X)){
      X[[i]] <- standardize_X(X[[i]])
    }
    
    # update data matrices by removing the latent variable - only for X and Y
    zi <- Z_std %*% zj
    xi_list <- lapply(seq_along(X), function(i) {
      X[[i]] %*% xj_list[[i]]
    })
    
    # deflation
    for(i in seq_along(X)) {
      
      xi <- xi_list[[i]]
      Xi <- X[[i]]
      
      reg <- apply(Xi, 2, function(x) {lm(x ~ xi)} )
      X[[i]] <- sapply(reg, function(x) {x[[2]]} )   
    }
    
    all.lambdas <- cbind(all.lambdas, lambdas)
    all.corr <- c(all.corr, corr.scca)
    for (i in seq_len(n_sets)){
      cc3.weight[[i]] <- cbind(cc3.weight[[i]], xj_list[[i]])
      cc3.CV[[i]] <- cbind(cc3.CV[[i]],xi_list[[i]])
    }
    cc3.weight.z <- cbind(cc3.weight.z, zj)
    cc3.CV.z <- cbind(cc3.CV.z, zi)
    
    canVar <- canVar + 1
  } # while canVar
  
  return(list(cc3.weight.x = cc3.weight,
              cc3.weight.z = cc3.weight.z,
              cc3.CV.x = cc3.CV,
              cc3.CV.z = cc3.CV.z,
              corr = all.corr,
              lambda = all.lambdas,
              num.CV = numCV))
}
