\name{getCCA3}
\alias{getCCA3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate supervised CCA on three data matrices}
\description{
%  ~~ A concise (1-5 lines) description of what the function does. ~~
getCCA3 starts the supervised sparse CCA with elastic net for 2 biological data sets and 1 design data set Z.
Data sets X and Y are strongly regularized (ridge regression); Z if neccessary (det(Z)==0).

}
\usage{
getCCA3(X,Y,Z, end=c(0.3,0.5,3), step=c(0.01,0.01,0.2), numCV=10,
n.r=10, max.counter.test=10)

}
\arguments{
  \item{X}{
    biological data set with dim(X): n rows, p features
}
  \item{Y}{
    biological data set with dim(Y): n rows, q features
}
  \item{Z}{
    design data set with dim(Z): n rows, r features; describes experimental design with binary vectors
}

  \item{end}{
    vector containing the three maximal lambda-values for the grid search
}
  \item{step}{
     vector containing the step size for grid search for the lambdas
}

  \item{numCV}{
    number of canonical variables, which will be calculated
  }

  \item{n.r}{
    number of resampling runs, 10 is okay
  }
  \item{max.counter.test}{
    number of random start vectors 
  }
}
\details{
%  ~~ If necessary, more details than the description above ~~
}
\value{
    \item{
        cc3.weight.x, cc3.weight.y, cc3.weight.z}{Matrix with weights of canonical variables for data set X, Y and Z respectively (columnwise)
    }
    \item{
        ccc3.CV.x, ccc3.CV.y, ccc3.CV.z }{Matrix with canonical variables for data set X, Y and Z respectively (columnwise)
    }
    \item{corr}{
        Vector with absolute correlation coefficients of canonical variables (columnwise): (corr(cv(X),cv(Z))+ corr(cv(Y),cv(Z))/2 
    }
    \item{lambda}{
        Matrix with 3 parameters lambda for each data set and each canonical variable (columns)
    }
    \item{numCV}{
        Number of canonical variables
    }
}
\author{
    Andrea Thum
}
\note{
%  ~~further notes~~
}

% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ supervised }
\keyword{ CCA }
\keyword{ sparse }% __ONLY ONE__ keyword per line
