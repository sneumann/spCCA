\name{save.CCA}
\alias{save.CCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate supervised CCA on three data matrices}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
   Store results of current canonical variable in a file. First row: information about lambdas and correlation coefficients. Next rows for all three data sets for each weight > 0 weight and feature name. Weights are sorted in decreasing order, so most important features are on top.
}
\usage{
save.CCA(CCA3,filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CCA3}{
        list with results from getCCA3
}
  \item{filename}{
        filename for saving results
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
}
\author{
    Andrea Thum
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

