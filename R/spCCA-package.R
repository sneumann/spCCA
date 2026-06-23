#' @keywords internal 
#' @title spCCA
#' @description This package provides functions to perform sparse canonical correlation analysis (spCCA) on two or more sets of variables and a design matrix as supervision.
#' @details The spCCA method is a sparse version of canonical correlation analysis (CCA) that can be used to identify linear combinations of variables from two or more sets that are maximally correlated. The spCCA method is particularly useful when the number of variables in each set is large and the number of observations is small. The spCCA method is implemented using the alternating direction method of multipliers (ADMM) algorithm. The package also provides functions to visualize the results of the spCCA analysis.
#' @name spCCA-package
#' @docType package
#' @keywords package
#' @version 0.1.7
#' @date 2026-06-03
#' @license GPL (>= 2)
#' @author 
#'   Anjana Bhat
#'   Andrea Thum
#'   Steffen Neumann
#'   Stefan Posch
#' @maintainer Steffen Neumann <sneumann@ipb-halle.de>
#' @references
#'   [@Hotel36]
#' @suggests testthat BiocStyle knitr rmarkdown 
#' @VignetteBuilder knitr 
"_PACKAGE"
