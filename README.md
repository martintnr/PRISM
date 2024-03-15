
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PleioVar

<!-- badges: start -->
<!-- badges: end -->

PleioVar takes as input GWAS data and outputs pleiotropic labels for a
list of variants.

## Installation

You can install the current version of PleioVar like so:

``` r
  if(!"devtools" %in% installed.packages()[, "Package"]) {
  install.packages("devtools") }

  if(!"PleioVar" %in% installed.packages()[, "Package"]) {
  devtools::install_github("martintnr/PleioVar") }
```

## Example

This is a basic example to obtain pleiotropic labels from GWAS summary
statistics.  
Parameters obtained from LHC-MR, LDscores of all variants, and simulated
GWAS ZScores are already included.  
You should specify the NbCores parameters if your computer can handle
parallel computations.

``` r
library(data.table)
library(dplyr)
library(stats)
library(stringr)
library(parallel)
library(matrixStats)
library(PleioVar)


if(!file.exists("PleioVar_example/")){system("mkdir PleioVar_example")}
if(!grepl("PleioVar_example", getwd(), fixed = TRUE)){setwd("PleioVar_example")}

Prepare_example_data(gzip = T)

ParametersTable <- fread("Data/ParametersTable.csv", header = T, sep = ",")
Index <- fread("Data/Index.csv", header = T, sep = ",")

ListofTraits <- unique(c(ParametersTable$X, ParametersTable$Y))

PleioVar_main(ListofTraits, ParametersTable, Index , NbCores = 1, gzip = F, keepIntermediateFiles = F)
```

In the Results/ folder can be found, for each trait, a file with
variants, p-values from PleioVar, and pleiotropy annotation.

``` r
library(ggplot2)

Graph <- Example_graph(ListofTraits, ParametersTable, Trait = "B4")

print(Graph)
```
