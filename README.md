# Low level infrastructure to handle MS spectra

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/sneumann/spCCA/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/sneumann/spCCA/actions?query=workflow%3AR-CMD-check-bioc)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

The spCCA package combines sparse (penalized) CCA for two biological data sets 
  with a third data set, the experimental design and tries to find correlations 
  between biological data sets and design data set, that are well-explainable.

For more information see the package
[homepage](https://sneumann.github.io/spCCA/).


# Installation

The package can be installed with

```r
install.packages("devtools")
devtools::install_github("sneumann/spCCA")
```

