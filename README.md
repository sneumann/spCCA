# supervised penalized CCA for n+1 data sets (n biological, 1 design)

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/sneumann/spCCA/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/sneumann/spCCA/actions?query=workflow%3AR-CMD-check-bioc)
[![GitHub License](https://img.shields.io/github/license/sneumann/spCCA)](https://opensource.org/license/gpl-2-0)


## Overview
This package implements supervised penalized canonical correlation analysis ( spCCA) to integrate multiple biological datasets/ multi-omics data with an experimental design encoded as a binary matrix.  
It aims to identify interpretable correlations between biological features and experimental conditions, making the relationships more biologically meaningful and easier to interpret.

For more information see the package
[homepage](https://sneumann.github.io/spCCA/).

---

## Features

- Sparse canonical correlation analysis (CCA)
- Integration of multiple biological datasets
- Incorporation of experimental design as a supervisory signal to guide the integration
- Improved interpretability through Elastic-Net penalization (feature selection)
- Enhanced computational efficiency through parallel processing

---

## Installation
The package can be installed with

```r
install.packages("devtools")
devtools::install_github("AnjuBhat247/spCCA")
```