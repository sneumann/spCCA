#' @keywords internal 
#' @title spCCA
#' @description This package provides functions to perform sparse canonical correlation analysis (spCCA) on two sets of variables.
#' @details The spCCA method is a sparse version of canonical correlation analysis (CCA) that can be used to identify linear combinations of variables from two sets that are maximally correlated. The spCCA method is particularly useful when the number of variables in each set is large and the number of observations is small. The spCCA method is implemented using the alternating direction method of multipliers (ADMM) algorithm. The package also provides functions to visualize the results of the spCCA analysis.
#' @name spCCA-package
#' @docType package
#' @keywords package
#' @version 0.1.16
#' @date 2024-12-25
#' @license GPL (>= 2)
#' @author 
#'   Andrea Thum
#'   Steffen Neumann
#'   Stefan Posch
#' @maintainer Steffen Neumann <sneumann@ipb-halle.de>
#' @references Please see https://arxiv.org/abs/1405.1534 for a manuscript describing spCCA.
#' @suggests testthat BiocStyle knitr rmarkdown BiocStyle
#' @VignetteBuilder knitr 
"_PACKAGE"
