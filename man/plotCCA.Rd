\name{plotCCA}
\alias{plotCCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate supervised CCA on three data matrices}
\description{
    Plot canonical variable and top two features for this variable for data sets X and Y
}
\usage{



plotCCA(CCA3, X, Y,filename,legend)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CCA3}{
        list with results from getCCA3
}
  \item{X, Y}{ 
        original biological data matrices without subtraction of canonical variables, for plots of top features
} 
  \item{filename}{CCA-project name, used to create filename}
  \item{legend}{names of data sets X, Y, Z for legend, e.g. c("Genes", "Metabolites","Design Patterns")}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
}
\author{
    Andrea Thum
}
\note{
    
}

\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

