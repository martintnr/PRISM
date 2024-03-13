
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PleioVar

<!-- badges: start -->
<!-- badges: end -->

PleioVar takes as input GWAS data and outputs pleiotropic labels for
HapMap3 variants.

## Installation

You can install the current version of PleioVar like so:

``` r
#install.packages("devtools")
devtools::install_github("martintnr/PleioVar")
library(PleioVar)
```

## Example

This is a basic example to obtain pleiotropic labels from simulated GWAS
data. Parameters obtained from LHC-MR are already included.

``` r

# If you want to process your own data with PleioVar, NewData needs to be True
NewData = F

 pkgs = c("data.table", "dplyr", "devtools","stats", "stringr", "parallel")
  pkgs.na = pkgs[!pkgs %in% installed.packages()[, "Package"]]
  
  if (length(pkgs.na) > 0) {
    install.packages(pkgs.na)
  }
  
  
  if(!"lhcMR" %in% installed.packages()[, "Package"] & NewData == T) {
  devtools::install_github("LizaDarrous/lhcMR")
  }

  if(!"PleioVar" %in% installed.packages()[, "Package"]) {
  devtools::install_github("martintnr/PleioVar")
  }


library(data.table)
library(dplyr)
library(stats)
library(stringr)
library(parallel)
library(lhcMR)
library(PleioVar)

setwd("/home/martin/Script/PleioVar/test")

Prepare_example_data()

PleioVar_main(ListofTraits, ParametersTable, NbCores = 1)
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
