
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PRISM: Pleiotropic Relationships to Infer the SNP Model

<!-- badges: start -->
<!-- badges: end -->

## Principle of PRISM

![PRISM](Github_Fig.png) PRISM takes as input GWAS summary statistics of
multiple traits and outputs significant labeled variant-trait effects.
Further information can be found in our
[pre-print](https://doi.org/10.1101/2024.06.01.24308193).

## Results visualization

You can see PRISM results and the causal network of any variant on our
[Shiny visualization tool](https://verbam01.shinyapps.io/PRISM/).

## Installation

To install the current version of PRISM, which should take a few
minutes:

``` r
  if(!"devtools" %in% installed.packages()[, "Package"]) {
  install.packages("devtools") }

  if(!"PRISM" %in% installed.packages()[, "Package"]) {
  devtools::install_github("martintnr/PRISM") }
```

PRISM has been tested on Linux: Ubuntu 22.04.4 LTS and macOS: Sonoma
14.5

## Example

This is a basic example to obtain labels from simulated GWAS summary
statistics.  
Parameters obtained from LHC-MR, LDscores of all variants, and simulated
GWAS standardized effect sizes are already included.  
You should specify the NbCores parameters if your computer can handle
parallel computations.

``` r
library(data.table)
library(dplyr)
library(stats)
library(stringr)
library(parallel)
library(matrixStats)
library(PRISM)


if(!grepl("PRISM_example", getwd(), fixed = TRUE)){
if(!file.exists("PRISM_example/")){system("mkdir PRISM_example")}
setwd("PRISM_example")}

Prepare_example_data(gzip = T)

ParametersTable <- fread("Data/ParametersTable.csv", header = T, sep = ",")
Index <- fread("Data/Index.csv", header = T, sep = ",")

ListofTraits <- unique(c(ParametersTable$X, ParametersTable$Y))

PRISM_main(ListofTraits, ParametersTable, Index , sourceGWAS = "Data/",
              NbCores = 1, gzip = F, keepIntermediateFiles = F)
```

On a single core, `PRISM_main` is expected to complete execution within
approximately 20 minutes. Upon completion, the Results/ directory will
contain, for each trait, a file that includes variants, p-values from
PRISM, and pleiotropy annotations.

To visualize a result [example
graph](https://github.com/martintnr/PRISM/blob/main/Example_graph_output.png):

``` r
library(ggplot2)

Graph <- Example_graph(ListofTraits, ParametersTable, Trait = "B4")

print(Graph)
```

To visualize an [example
network](https://github.com/martintnr/PRISM/blob/main/Example_network_output.png),
with Shiny:

``` r
library(splitstackshape)
library(visNetwork)
library(shiny)

Network <- Example_network(ListofTraits, ParametersTable, Variant = "1:5341323:G:A")
```

If you want to use PRISM on real data, you can use this
[vignette](https://github.com/martintnr/PRISM/blob/main/vignettes/PRISM_vignette.md).
