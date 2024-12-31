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
#' @param X,Y,Z three data sets, Z is design data set, (STN: all three ?) already normalized to zero mean and variance 1
#' @param start,end,step start, end and stepsize for lambda.{x,y,z} grid search.
#' @param n.r number of resampling steps
#' @param max.counter.test number of random start vectors for iteration (inner loop)
#'
#' @return list with three components:
#'       best.lambdas: optimal lambdas for the three data sets
#        corr: best test set correlation for best vector
#        bestVector: all three eigenvectors for x, y, z with best correlation concatenated
#' 
#' @export
#' @importFrom stats var cor medien runif lm
#' @importFrom graphics axis lines
#' @importFrom grDevices dev.off pdf rainbow 
#' 
#' @author Andrea Thum, Elena Parkhomenko
#' @examples
#' TRUE
get.best.lambdas <- function(X, Y, Z,
                             start = c(0, 0, 0),
                             end = c(0.2, 0.2, 2),
                             step = c(0.02, 0.02, 0.2),
                             n.r = 10,
                             max.counter.test = 10)
{
  lambda.x.seq <- seq(start[1], end[1], by = step[1])
  lambda.y.seq <- seq(start[2], end[2], by = step[2])
  lambda.z.seq <- seq(start[3], end[3], by = step[3])
  
  n.lambdas.x <-  length(lambda.x.seq)
  n.lambdas.y <-  length(lambda.y.seq)
  n.lambdas.z <-  length(lambda.z.seq)
  
  p <- dim(X)[2]
  q <- dim(Y)[2]
  r <- dim(Z)[2]
  
  dims <- c(p, q, r)
  n.sample <- dim(X)[1]
  n.r.sample <- max(trunc(n.sample / 8), 3)#at least 1/8 of the features - otherwise problems with Z may occur in testing sample
  whole.sample <- seq(1, n.sample)
  
  Xp <- (diag(1 / sqrt(diag(var(X))))) %*% t(X)#strong regularization
  Yp <- (diag(1 / sqrt(diag(var(Y))))) %*% t(Y)
  if (abs(det(var(Z))) > 10 ^ -20) {
    #Zp <- solve(var(Z)) %*% t(Z)
    Zp <- ginv(var(Z)) %*% t(Z)
  } else {
    # if var(Z) is invertible
    Zp <- (diag(1 / sqrt(diag(var(Z))))) %*% t(Z) # otherwise: regularize
  }
  
  # Pseudomatrices * ...
  XpZ <- Xp %*% Z
  YpZ <- Yp %*% Z
  ZpX <- Zp %*% X
  ZpY <- Zp %*% Y
  
  maxIteration <- 500
  
  test.corr.scca <- array(0, dim = c(n.lambdas.x, n.lambdas.y, n.lambdas.z)) # corr. for test data set with trained vectors
  canVarVector <-  array(0, dim = c(n.lambdas.x, n.lambdas.y, n.lambdas.z, p + q + r))
  
  # while-loop here, because test data sets can lead to 0-variation in special cases
  # lists of pseudo inverse matrices
  x.train.list <- list()
  y.train.list <- list()
  z.train.list <- list()
  
  x.test.list <- list()
  y.test.list <- list()
  z.test.list <- list()
  
  xpz.predict.list <- list()
  ypz.predict.list <- list()
  zpx.predict.list <- list()
  zpy.predict.list <- list()
  
  # produce test and training data sets
  k = 1
  while (k <= n.r) {
    testing.sample = c()
    training.sample <- whole.sample
    
    for (i in 1:n.r.sample) {
      x = training.sample[round(runif(1, 1, length(training.sample)))]
      training.sample = training.sample[training.sample != x]
      testing.sample <- c(testing.sample, x)
    }
    
    x.predict <- X[training.sample, ]
    y.predict <- Y[training.sample, ]
    z.predict <- Z[training.sample, ]
    
    x.test <- X[testing.sample, ]
    y.test <- Y[testing.sample, ]
    z.test <- Z[testing.sample, ]
    
    if (max(abs(var(z.test))) == 0)
      next # only zeros, try again
    
    # Standardize data
    x.predict <- apply(x.predict, 2, function(y) {
      if (var(y) == 0)
        y - mean(y)
      else
        (y - mean(y)) / var(y)
    })
    y.predict <- apply(y.predict, 2, function(y) {
      if (var(y) == 0)
        y - mean(y)
      else
        (y - mean(y)) / var(y)
    })
    z.predict <- apply(z.predict, 2, function(y) {
      if (var(y) == 0)
        y - mean(y)
      else
        (y - mean(y)) / var(y)
    })
    
    x.test <- apply(x.test, 2, function(y) {
      if (var(y) == 0)
        y - mean(y)
      else
        (y - mean(y)) / var(y)
    })
    y.test <- apply(y.test, 2, function(y) {
      if (var(y) == 0)
        y - mean(y)
      else
        (y - mean(y)) / var(y)
    })
    z.test <- apply(z.test, 2, function(y) {
      if (var(y) == 0)
        y - mean(y)
      else
        (y - mean(y)) / var(y)
    })
    
    # Pseudoinverse
    Xp.pred <- (diag(1 / sqrt(diag(var(x.predict ))))) %*% t(x.predict) # strong regularization
    Yp.pred <- (diag(1 / sqrt(diag(var(y.predict ))))) %*% t(y.predict)
    
    if (abs(det(var(z.predict))) > 10 ^ -20) {
      #Zp.pred <- solve(var(z.predict)) %*% t(z.predict)  #if var(Z) is invertible
      Zp.pred <- ginv(var(z.predict)) %*% t(z.predict)  #if var(Z) is invertible
    } else{
      Zp.pred <- (diag(1 / sqrt(diag(var(z.predict))))) %*% t(z.predict) # otherwise: regularize
    }
    
    xpz.predict = Xp.pred %*% z.predict
    ypz.predict = Yp.pred %*% z.predict
    zpx.predict = Zp.pred %*% x.predict
    zpy.predict = Zp.pred %*% y.predict
    
    x.test.list[[k]] <- x.test
    y.test.list[[k]] <- y.test
    z.test.list[[k]] <- z.test
    
    x.train.list[[k]] <- x.predict
    y.train.list[[k]] <- y.predict
    z.train.list[[k]] <- z.predict
    
    xpz.predict.list[[k]] <- xpz.predict
    ypz.predict.list[[k]] <- ypz.predict
    zpx.predict.list[[k]] <- zpx.predict
    zpy.predict.list[[k]] <- zpy.predict
    
    k <- k + 1
    
  } # while k
  
  #now cross-val
  for (j.lambda.y in 1:n.lambdas.y)
  {
    lambda.y <- lambda.y.seq[j.lambda.y]
    for (j.lambda.x in 1:n.lambdas.x)
    {
      lambda.x <- lambda.x.seq[j.lambda.x]
      for (j.lambda.z in 1:n.lambdas.z)
      {
        lambda.z <- lambda.z.seq[j.lambda.z]
        if (is.na(test.corr.scca[j.lambda.x, j.lambda.y, j.lambda.z]))
          next

        all.Patterns <- list()
        i.r = 0
        for (i.r in 1:n.r) {
          # Resampling
          best.corr.train <- NA
          for (counter.test in 1:max.counter.test) {
            uv <- scca.function3Z(XpZ=xpz.predict.list[[i.r]],
                                  YpZ=ypz.predict.list[[i.r]],
                                  ZpX=zpx.predict.list[[i.r]],
                                  ZpY=zpy.predict.list[[i.r]],
                                  dims=dims,
                                  lambda.x=lambda.x, lambda.y=lambda.y, lambda.z=lambda.z,
                                  max.iter=maxIteration,
                                  shift=2) # z.initial=bestComb[[4]])
            
            if (uv$null ||
                uv$i == maxIteration)
              next
            # lambda too big for prediction data set
            if (is.na(best.corr.train))
              best.corr.train = 0 # converged for first time
            
            xj <- uv$x.new	# sparse singular vector (canonical vector for Y)
            yj <- uv$y.new	# sparse singular vector (canonical vector for X)
            zj <- uv$z.new
            
            if (var(x.test.list[[i.r]] %*% xj) == 0 ||
                var(y.test.list[[i.r]] %*% yj) == 0 ||
                var(z.test.list[[i.r]] %*% zj) == 0)
              next
            else
              corr.train <- sum(c(abs(
                cor(x.train.list[[i.r]] %*% xj, z.train.list[[i.r]] %*% zj)), 
                abs(cor(y.train.list[[i.r]] %*% yj, z.train.list[[i.r]] %*% zj) ))) / 2
            if (corr.train > best.corr.train) {
              best.corr.train <- corr.train
              xyz.vector <- c(xj, yj, zj)
            }
            
          } # for counter.test
          if (is.na(best.corr.train)) {
            # lambda too big
            test.corr.scca[j.lambda.x, j.lambda.y, j.lambda.z] <- NA
            all.Patterns <- list()
            
            break
          }
          
          test.corr.scca[j.lambda.x, j.lambda.y, j.lambda.z] <-
            test.corr.scca[j.lambda.x, j.lambda.y, j.lambda.z] + best.corr.train
          
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
          # Following two lines commented out by athum: 
          # all.Patterns<-addPattern(all.Patterns,best.corr.train,xyz.vector,xyz.vector[(p+q+1):(p+q+r)])
          # Vectors<-cbind(Vectors,xyz.vector)
        } # while resampling
        
        # cluster xyz.vectors
        if (length(all.Patterns) == 0) {
          test.corr.scca[j.lambda.x, j.lambda.y, j.lambda.z] <- NA
        } else {
          ResamplingCorr <- unlist(lapply(all.Patterns, 
                                          function(item) {item[[1]]} ))
          
          # determine median
          ResamplingCorr2 <- ResamplingCorr
          if (n.r %% 2 == 0)
            ResamplingCorr2 <- sort(ResamplingCorr)[2:length(ResamplingCorr)]
          
          wRC <- which(ResamplingCorr == median(ResamplingCorr2))[1]

          canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, ] <- rowMeans(all.Patterns[[wRC]][[2]])
          
          cx <- canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, 1:dims[1]]
          canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, 1:dims[1]] <-
            canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, 1:dims[1]] / (sqrt(t(cx) %*% cx))
          
          cy <- canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, 
                             (dims[1] + 1):(dims[1] + dims[2])]
          canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, 
                       (dims[1] + 1):(dims[1] + dims[2])] <- 
            canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, 
                         (dims[1] + 1):(dims[1] + dims[2])] / (sqrt(t(cy) %*% cy))
          
          cz <- canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, 
                             (dims[1] + dims[2] + 1):(dims[1] + dims[2] + dims[3])]
          canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, 
                       (dims[1] + dims[2] + 1):(dims[1] + dims[2] + dims[3])] <- 
            canVarVector[j.lambda.x, j.lambda.y, j.lambda.z, 
                         (dims[1] + dims[2] + 1):(dims[1] + dims[2] + dims[3])] / (sqrt(t(cz) %*% cz))
        }
        
      } # lambda.z
    } # lambda.x
  } # lambda.y
  
  test.corr.scca[is.na(test.corr.scca)] <- 0
  test.corr.scca <- test.corr.scca / n.r
  max.corr <- max(test.corr.scca)
  if (max.corr == 0)
    return(NULL) # no correlations found
  
  i.lambda <- which(test.corr.scca == max.corr, TRUE)[1, ]
  
  lambda.x = i.lambda[1] * step[1] - step[1]
  lambda.y = i.lambda[2] * step[2] - step[2]
  lambda.z = i.lambda[3] * step[3] - step[3]
  bestVector <- canVarVector[i.lambda[1], i.lambda[2], i.lambda[3], ]
  
  return(list(best.lambda.x = lambda.x,
              best.lambda.y = lambda.y,
              best.lambda.z = lambda.z,
              corr = max.corr,
              bestVector = bestVector))
}
