\name{scca.function3Z}
\alias{scca.function3Z}
\title{Calculate supervised CCA on three data matrices}
\description{
    Implements the modified iterative power method for three (regularized) data sets with lasso penalty for given lasso parameters lambda.\{x,y,z\}.
    Pseudomatrices are used as input for faster calculation.
}
\usage{
scca.function3Z(XpZ,YpZ,ZpX,ZpY, dims,lambda.x, lambda.y, lambda.z,max.iter = 500, 
shift=2, x.initial=NULL, y.initial=NULL, z.initial=NULL)
}

\arguments{
  \item{XpZ, YpZ, ZpX, ZpY}{
        For regularized pseudomatrix Xp of X: regularize(t(X)*X)) * t(X); 
        then XpZ = Xp * Z is input. Z is the design data set.
      }
  \item{dims}{
        dims[1:3]: number of features of data sets X, Y, Z
      }
  \item{lambda.x, lambda.y, lambda.z}{
        lasso parameters lambda for data sets. lambda = 0: no sparsity.
      }
  \item{max.iter}{
        number of maximal iterations for power method. 
      }
  \item{shift}{
        start lasso regularization after shift iteration steps.
      }
  \item{x.initial, y.initial, z.initial}{
        initial vectors for iteration. If NULL: filled with random values
      }
    }
    
\value{
    \item{x.new}{weight vector for data set X}
    \item{y.new}{weight vector for data set Y}
    \item{z.new}{weight vector for data set Z}
    \item{i}{number of iterations steps until convergence}
    \item{null}{true if at least one lambda is too high and all weights of one weight vector are zero}
}

\author{
    Andrea Thum, Elena Parkhomenko
}
\note{
    Test output variable i, if iteration did converge (i<max.iter)!
}
\examples{

}


\keyword{ CCA }
\keyword{ lasso }
