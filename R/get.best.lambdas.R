#' Calculate supervised CCA on three data matrices
#' 
#' Calculate correlation coefficients for sparse CCA for all combinations 
#' of lambdas (using grid search).
#' 
#' Adequate numbers for the grid search depend on the data set size. 
#' As a starting point 0.0 is recommended, this is equivalent to no sparsity. 
#' Depending on the number of features, 
#' for about 1000 features, as ending point 0.25 should be high enough, 
#' for about 10 features it is more like 2. 
#' For the step size 0.01 is a good number for sets with 1000 features, 0.1 for 10 features. 
#' For a good resampling quality set n.r to about ten.
#' The higher the number of random start vectors, the better. 
#' Depending on the instability of the data set between 5 and 20 should be good.
#' 
#' As design data set Z is small, lambdas are higher and bigger steps can be used 
#' for data set with about 1000 features, lambda.end=0.3 should be high enough; 
#' use higher end values for smaller data sets and lower for bigger data sets.
#'
#' @param X,Z biological data sets, Z is design data set
#' @param lambda.grid set of lambda combinations
#' @param n.r number of resampling steps
#' @param max.counter.test number of random start vectors for iteration (inner loop)
#'
#' @return list with three components:
#'       best.lambdas: optimal lambdas for the data sets
#        corr: best test set correlation for best vector
#        bestVector: all eigenvectors for X and Z with best correlation concatenated
#' 
#' @export
#' @importFrom stats var cor median runif lm sd
#' @importFrom graphics axis lines
#' @importFrom grDevices dev.off pdf rainbow 
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG %dorng% registerDoRNG
#' 
#' @author Anjana Bhat, Andrea Thum, Elena Parkhomenko
#' @examples
#' TRUE
get.best.lambdas <- function(X, Z,
                             lambda.grid,
                             n.r = 10,
                             max.counter.test = 10)
{
  n.sets <- length(X)
  
  dims.X <- sapply(X, ncol)
  dims <- c(dims.X, ncol(Z))
  
  n.sample <- nrow(X[[1]])
  n.r.sample <- max(trunc(n.sample / 8), 3)#at least 1/8 of the features - otherwise problems with Z may occur in testing sample
  whole.sample <- seq(1, n.sample)
  
  maxIteration <- 100
  test.corr <- 0

  # while-loop here, because test data sets can lead to 0-variation in special cases
  # lists of pseudo inverse matrices
  X.train.list <- vector("list", n.r)
  X.test.list  <- vector("list", n.r)
  Z.train.list <- list()
  Z.test.list <- list()
  XpZ.train.list <- vector("list", n.r)
  ZpX.train.list <- vector("list", n.r)
  
  # produce test and training data sets
  k = 1
  while (k <= n.r) {
    
    testing.sample <- sample(whole.sample, n.r.sample)
    training.sample <- setdiff(whole.sample, testing.sample)
    
    X.train <- lapply(X, function(mat) mat[training.sample, , drop = FALSE])
    X.test  <- lapply(X, function(mat) mat[testing.sample, , drop = FALSE])
    
    Z.train <- Z[training.sample, , drop = FALSE]
    Z.test  <- Z[testing.sample, , drop = FALSE]
    
    
    if (max(abs(var(Z.test))) == 0)
      next # only zeros, try again
    
    # Standardize data: replaced var with sd. Also, independent standardization of test 
    # might lead to data leakage, hence we use mu and sdv of train set
    
    standardize.train <- function(M) {
      mu <- colMeans(M)
      sdv <- apply(M, 2, sd)
      sdv[sdv == 0] <- 1
      
      M.scaled <- sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
      list(X = M.scaled, mu = mu, sd = sdv)
    }
    
    standardize.test <- function(M, mu, sdv) {
      sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
    }
    
    
    for (i in seq_along(X.train)) {
      tmp <- standardize.train(X.train[[i]])
      X.train[[i]] <- tmp$X
      X.test[[i]]  <- standardize.test(X.test[[i]], tmp$mu, tmp$sd)
    }
    
    # design matrix Z might be full rank or singular. If columns are singular, we use Moore-Penrose pseudoinverse. If not, we compute inverse of the matrix
    Z.std <- apply(Z.train, 2, function(y){
      if (var(y) == 0)
        y - mean(y)
      else
        (y - mean(y)) / var(y)
    })
    
    # determine pseudo matrix once - Z does not change
    if (abs(det(var(Z.std))) > 10 ^ -20) {
      # if var(Z) is invertible
      Zp.train <- ginv(var(Z.std)) %*% t(Z.std) # cov(Z) proportional to Z.T *Z
      
    } else{
      Zp.train <- (diag(1 / sqrt(diag(var(Z.std))))) %*% t(Z.std) # otherwise: regularize
    }
    
    # standardizing Z-train now
    tmp <- standardize.train(Z.train)
    Z.train <- tmp$X
    Z.test  <- standardize.test(Z.test, tmp$mu, tmp$sd)
    
    # Pseudomatrix 
    Xp.train <- lapply(X.train, t)
    
    XpZ.train <- list()
    ZpX.train <- list()
    
    for (i in seq_along(Xp.train)){
      XpZ.train[[i]] <- Xp.train[[i]] %*% Z.train
      ZpX.train[[i]] <- Zp.train %*% X.train[[i]]
    }
    
    X.test.list[[k]] <- X.test
    Z.test.list[[k]] <- Z.test
    
    X.train.list[[k]] <- X.train
    Z.train.list[[k]] <- Z.train
    
    XpZ.train.list[[k]] <- XpZ.train
    ZpX.train.list[[k]] <- ZpX.train
    
    k <- k + 1
    
  } # while k
  
  # now cross-validate: parallel processing
  n_cores <- detectCores()
  workers <- max(2, n_cores-2)
  # Limit to 2 when CRAN is checking
  if (identical(tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", "false")), "true")) {
    workers <- 2
  }
  cl <- makeCluster(workers)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    options(error = function() {
      traceback(2)
      stop("Worker error")
    })
  })
  registerDoRNG(3)
  
  
  results <- foreach(j = 1:nrow(lambda.grid),
                     .combine = 'c',
                     .export = "scca.function3Z",
                     .packages = c("MASS")) %dopar% {

               lambda.x <- as.numeric(lambda.grid[j, 1:n.sets]) 
               lambda.z <- as.numeric(lambda.grid[j, n.sets+1])
               all.Patterns <- list()
               test.corr <- 0 
               valid.count <- 0

              i.r = 0
              for (i.r in 1:n.r) {
                # Resampling
                best.corr.train <- NA
                prev.best <- -Inf
                stall <- 0
                stall.max <- 2
                tol <- 1e-4
                
                for (counter.test in 1:max.counter.test) {
                  uv <- scca.function3Z(XpZ=XpZ.train.list[[i.r]],
                                        ZpX=ZpX.train.list[[i.r]],
                                        dims=dims,
                                        lambda.x=lambda.x, lambda.z=lambda.z,
                                        max.iter=maxIteration,
                                        shift=2) # z.initial=bestComb[[4]])
                  
                  if (uv$null || uv$n== maxIteration)
                    next
                  # lambda too big for prediction data set
                  if (is.na(best.corr.train))
                    best.corr.train = 0 # converged for first time
                  
                  xj <- uv$x.new	# sparse singular vector (canonical vector for X)
                  zj <- uv$z.new
                  
                  if (any(mapply(function(x_t, x_j) var(x_t %*% x_j) == 0, X.test.list[[i.r]], xj )) ||
                      var(Z.test.list[[i.r]] %*% zj) == 0)
                    next
                  else
                    corr.train <- mean(mapply(function(dat, vec) {
                      abs(cor(
                        dat %*% vec, Z.train.list[[i.r]] %*% zj
                      ))
                    }, X.train.list[[i.r]], xj))
                  
                  if (corr.train > best.corr.train) {
                    best.corr.train <- corr.train
                    xyz.vector <- unlist(c(xj, list(zj)))
                  }
                  
                  if(best.corr.train <= prev.best + tol) {
                    stall <- stall + 1
                  } else {
                    stall <- 0
                    prev.best <- best.corr.train
                  }
                  
                  if(stall >= stall.max) {
                    break
                  }
                } # for counter.test
                if (is.na(best.corr.train)) {
                  # lambda too big
                  test.corr <- NA
                  all.Patterns <- list()
                  valid.count <- 0
                  break
                }
                
                test.corr <-  test.corr + best.corr.train
                valid.count <- valid.count + 1
                
                if (length(all.Patterns) == 0) {
                  all.Patterns <- list(list(best.corr.train,
                                            matrix(xyz.vector,
                                                   ncol = 1,
                                                   nrow = length(xyz.vector))))
                } else {
                  all.Patterns <- append(all.Patterns, 
                                         list(list(best.corr.train,
                                                   matrix(xyz.vector,
                                                          ncol = 1,
                                                          nrow = length(xyz.vector)))))
                }
                
              } # while resampling
              
              # cluster xyz.vectors
              if (length(all.Patterns) == 0) {
                test.corr <- NA
                canVarVector <- NA
              } else {
                ResamplingCorr <- unlist(lapply(all.Patterns, 
                                                function(item) {item[[1]]} ))
                
                # determine median
                ResamplingCorr2 <- ResamplingCorr
                if (n.r %% 2 == 0)
                  ResamplingCorr2 <- sort(ResamplingCorr)[2:length(ResamplingCorr)]
                
                wRC <- which(ResamplingCorr == median(ResamplingCorr2))[1]
      
                canVarVector <- rowMeans(all.Patterns[[wRC]][[2]])
                
                starts <- c(1, cumsum(dims)[-length(dims)] + 1)
                ends <- cumsum(dims)
                
                for(k in seq_along(dims)) {
                  
                  ck <- canVarVector[starts[k]:ends[k]]
                  canVarVector[starts[k]:ends[k]] <- ck / c(sqrt(t(ck) %*% ck))
                }   
                test.corr <- test.corr / valid.count
              } 
              list(list(lambda = lambda.grid[j, ],
                        test.corr = test.corr,
                        canVarVector = canVarVector
              ))
              
  } # foreach end
  on.exit(stopCluster(cl))
  test.corr.scca <- sapply(results, `[[`, "test.corr")
  if(all(is.na(test.corr.scca)))
    return(NULL) # no correlations found
  max.corr <- max(test.corr.scca, na.rm = TRUE)
  i.lambda <- which.max(test.corr.scca)
    
  
  best.lambda <- results[[i.lambda]][[1]]
  lambda.x = best.lambda[1:n.sets]
  lambda.z = best.lambda[n.sets+1]
  bestVector <- results[[i.lambda]][[3]]
  
  return(list(best.lambda.x = lambda.x,
              best.lambda.z = lambda.z,
              corr = max.corr,
              bestVector = bestVector))
}
