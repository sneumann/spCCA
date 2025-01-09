\name{spCCA-package}
\alias{spCCA-package}
\alias{spCCA}
\docType{package}
\encoding{utf8}     
\title{
Supervised Sparse Canonical Correlation Analysis
}
\description{
  Obtain Canonical Variables that explain correlations between data
  matrices.  
}
\details{
\tabular{ll}{
Package: \tab spCCA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-06-15\cr
License: \tab under GPL license\cr
}
Also check the Vignette
}
\author{
  Andrea Thum, Steffen Neumann, Stefan Posch

Maintainer: andrea.thum@informatik.uni-halle.de
}
\references{

1. Supervised Penalized Canonical Correlation Analysis, Andrea Thum, Lore Westphal, Tilo Lübken, Sabine Rosahl, Steffen Neumann and Stefan Posch, unpublished

2. A learning algorithm for adaptive canonical correlation analysis of several data sets, Javier Vía and Ignacio Santamaría and Jesús Pérez,
Neural Networks 2006, volume 20,  pp. 139-152

3. Correlating multiple SNPs and multiple disease phenotypes: penalized non-linear canonical correlation analysis, Sandra Waaijenborg and Aeilko H. Zwinderman, Bioinformatics 2009, volume 25 (21),  pp. 2764-2771

4. Relations between two sets of variates, Harold Hotelling, Biometrika 1936, volume 28, pp. 321-377



}
\keyword{ package, CCA, biostatistics, sparse, supervised }
\seealso{
}
\examples{

X=as.matrix(read.delim(system.file("Example/simX.txt", package = "spCCA")))
Y=as.matrix(read.delim(system.file("Example/simY.txt", package = "spCCA")))
Z=as.matrix(read.delim(system.file("Example/simZ.txt", package = "spCCA")))

CCA3 <- getCCA3(X,Y,Z,numCV=2,  max.counter.test=1, n.r=1)


}
