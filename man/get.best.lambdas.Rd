\name{get.best.lambdas}
\alias{get.best.lambdas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate supervised CCA on three data matrices}
\description{
    Calculate correlaion coefficients for sparse CCA for all combinations of lambdas (grid search). 
}
\usage{
get.best.lambdas(X,Y,Z, start=c(0,0,0), end=c(0.2,0.2,2), step=c(0.02,0.02,0.2), 
n.r=10, max.counter.test=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X, Y, Z}{
        the data matrices, Z is supposed to be the design data set
}
  \item{start}{
        vector with three numbers for data sets X, Y and Z; starting points for grid search for lambdas.
}
  \item{end}{
        vector with three numbers for data sets X, Y and Z; ending points for grid search for lambdas. 
}
  \item{step}{
        vector with three numbers for data sets X, Y and Z; step size for grid search for lambdas
}
  \item{n.r}{
        number of resampling runs
}   
  \item{max.counter.test}{
        number of random start vectors for iteration
}

}
\details{
    Adequate numbers for the grid search depend on the data set size. As a starting point 0.0 is recommended, this is equivalent to no sparsity. Depending on the number of features, for about 1000 features, as ending point 0.25 should be high enough, for about 10 features it is more like 2. For the step size 0.01 is a good number for sets with 1000 features, 0.1 for 10 features. 

For a good resampling quality set n.r to about ten.

The higher the number of random start vectors, the better. Depending on the instability of the data set between 5 and 20 should be good.
}
\value{
    \item{best.lambda.x, best.lambda.y, best.lambda.z}{
       optimal lambdas for the three data sets.
    }
    \item{corr}{
        best test set correlation for best vector
    } 
    \item{bestVector}{
        vector; all three eigenvectors for x, y, z with best correlation concatenated
    }
}
\author{
    Andrea Thum, Elena Parkhomenko
}
\note{
    
}

\examples{

}

\seealso{
}



