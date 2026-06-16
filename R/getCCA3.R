#' Calculate supervised CCA on three or more data matrices
#'
#' getCCA3() starts the supervised sparse CCA 
#' 
#' 
#' 
#' 
#' The getCCA3() function performs sparse 
#' canonical correlation analysis on three or more data matrices 
#' (n Biological data sets in X, and one design data set Z) with elastic net. 
#' Data sets in X are strongly regularized (ridge regression); Z if neccessary (det(Z)==0).
#' 
#' The getCCA3 first normalizes the data matrices (column-wise mean centering and unit variance), 
#' computes canonical variables, 
#' and updates the data matrices by removing the latent variables iteratively 
#' for a specified number of canonical variables (numCV). 
#' It returns a list containing matrices with weight vectors and canonical variables, 
#' correlation coefficients, and optimal sparsity parameters (lambda values) 
#' for each canonical variable.
#' 
#' @param X list of biological data sets with dim(X1, X2,..): n rows/samples, p features
#' @param Z design data set with dim(Z): n rows/samples, r features; describes experimental design with binary vectors
#' @param end,step end and stepsize for lambda.(x,y,z) grid search.
#' @param n.r number of resampling runs, 10 is okay
#' @param max.counter.test number of random start vectors for iteration (inner loop)
#' @param numCV number of canonical variables which will be calculated
#' @param grid.search search strategy for lambda: exhaustive considers all possible lambda combinations, random considers sampled set of lambda combinations
#' @param n.comb number of lambda combinations to be sampled if random grid search is selected
#'
#' @return A list with the elements 
#'   cc3.weights.(xz): Matrices with weights of canonical variables for data sets X and Z respectively (columnwise)
#'   cc3.CV.(xz): matrices with canonical variables for data set X and Z respectively (columnwise)
#'   corr: matrix with correlation coefficient for pairwise (X-Z,Y-Z) and correlation of can. variable: mean(abs(corr(X-Z)),abs(corr(Y-Z))) (columnwise)
#'   lambda: best lambda.(xyz) for each data set and for each can. var. (columnwise)
#'   numCV: Number of canonical variables
#'   
#' @keywords supervised CCA sparse
#' @author Anjana Bhat, Andrea Thum
#' @importFrom MASS ginv
#' @importFrom stats sd
#' @export
#'
#' @examples
#' TRUE
getCCA3 <- function(X, Z,
                    end, # c(0.3, 0.5, 3),
                    step, # c(0.01, 0.01, 0.2),
                    numCV = 10,
                    n.r = 10,
                    max.counter.test = 10,
                    grid.search = 'exhaustive',
                    n.comb = 5000) {
  if (length(end) != length(X)+1 || length(step) != length(X)+1)
    stop('Length of end or stepsize does not match the total number of datasets')
  
  n.values <- mapply(function(e, s) {
    
    floor(e / s) + 1
    
  }, end, step)
  total.combinations <- prod(n.values)
  
  if (total.combinations > 20000 && grid.search == "exhaustive")
    message("\033[1;31m Total lambda combinations = ", total.combinations,
    ". The current lambda grid is very large and may take a long time to evaluate. ",
    "For faster execution, consider stopping the run and using random search (default n.comb = 5000) or adjusting end - step values.\033[0m")
  
  n.sets <- length(X)
  lambda.x.seq <- vector("list",n.sets)
  for (i in 1:n.sets) {
    lambda.x.seq[[i]] <- seq(0, end[i], by = step[i])
  }
  lambda.z.seq <- seq(0, end[n.sets+1], by = step[n.sets+1])
  lambda.list <- c(lambda.x.seq, list(lambda.z.seq))
  lambda.grid.full <- expand.grid(lambda.list)
  
  if (grid.search == "exhaustive")
    lambda.grid <- lambda.grid.full
  else if (grid.search == "random")
    lambda.grid <- lambda.grid.full[
      sample(nrow(lambda.grid.full), n.comb),
    ]
  else
    stop("Argument 'grid.search' must be either 'exhaustive' or 'random'.")
  
  Total.data <- c(X,list(Z))
  
  if (length(unique(sapply(Total.data,nrow)))!=1)
    stop('Data Matrices have different sample sizes: ',sapply(Total.data,nrow))

  cc3.weight <- vector("list", n.sets)
  cc3.CV <- vector("list", n.sets)
  
  for (i in seq_len(n.sets)) { 
    cc3.weight[[i]] <- matrix(ncol = 0, nrow = ncol(X[[i]]))
    rownames(cc3.weight[[i]]) <- colnames(X[[i]])
    
    cc3.CV[[i]] <- matrix(ncol = 0, nrow = nrow(X[[i]]))
    rownames(cc3.CV[[i]]) <- rownames(X[[i]])
  }
  
  cc3.weight.z <- matrix(ncol = 0, nrow = dim(Z)[2])
  cc3.CV.z <- matrix(ncol = 0, nrow = dim(Z)[1])
  rownames(cc3.weight.z) <- colnames(Z)
  rownames(cc3.CV.z) <- rownames(Z)
  
  all.lambdas <- matrix(ncol = 0, nrow = length(Total.data))
  all.corr <- c()
  
  Z.std <- apply(Z, 2, function(y) {
    if (sd(y) == 0)
      y - mean(y)
    else
      (y - mean(y)) / sd(y)
  })
  
  dims <- sapply(Total.data, ncol)
  
  canVar = 1
  while (canVar <= numCV) {
   
    # get best combination of sparsity parameters
    results <- get.best.lambdas(X, Z,
      lambda.grid = lambda.grid,
      n.r = n.r,
      max.counter.test = max.counter.test) 
    
    if (is.null(results)) {
      break
    }
    lambdas <- c(results$best.lambda.x,
                 results$best.lambda.z)
    corr.scca <- results$corr
    
    # split vectors
    split.idx <- rep(seq_along(dims), dims)
    vectors <- split(results$bestVector, split.idx)
    
    zj <- vectors[[length(vectors)]]
    
    xj.list <- vectors[-length(vectors)]
    
    # standardize X before deflation
    standardize.X <- function(M){
      mu <- colMeans(M)
      sdv <- apply(M, 2, sd)
      sdv[sdv == 0] <- 1
      
      M.scaled <- sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
      M.scaled
    }
    for (i in seq_along(X)){
      X[[i]] <- standardize.X(X[[i]])
    }
    
    # update data matrices by removing the latent variable - only for X and Y
    zi <- Z.std %*% zj
    xi.list <- lapply(seq_along(X), function(i) {
      X[[i]] %*% xj.list[[i]]
    })
    
    # deflation
    for(i in seq_along(X)) {
      
      xi <- xi.list[[i]]
      Xi <- X[[i]]
      
      reg <- apply(Xi, 2, function(x) {lm(x ~ xi)} )
      X[[i]] <- sapply(reg, function(x) {x[[2]]} )   
    }
    
    all.lambdas <- cbind(all.lambdas, lambdas)
    all.corr <- c(all.corr, corr.scca)
    for (i in seq_len(n.sets)){
      cc3.weight[[i]] <- cbind(cc3.weight[[i]], xj.list[[i]])
      cc3.CV[[i]] <- cbind(cc3.CV[[i]],xi.list[[i]])
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
