
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PRISM: Pleiotropic Relationships to Infer the SNP Model

<!-- badges: start -->
<!-- badges: end -->

## Principle of PRISM

![PRISM](Github_Fig.png) PRISM takes as input GWAS summary statistics of
multiple traits and outputs significant labeled variant-trait effects, and a causal network for each significant genetic variant.

## Results visualization

You can see PRISM results and the causal network of any significant variant on our
[Shiny visualization tool](https://verbam01.shinyapps.io/PRISM/)

## Preprint

You can find our preprint: [Inferring genetic variant causal network by leveraging pleiotropy](https://www.medrxiv.org/content/10.1101/2024.06.01.24308193v1)

## Installation

You can install the current version of PRISM like so:

``` r
  if(!"devtools" %in% installed.packages()[, "Package"]) {
  install.packages("devtools") }

  if(!"PRISM" %in% installed.packages()[, "Package"]) {
  devtools::install_github("martintnr/PRISM") }
```

## Example

This is a basic example to obtain labels from simulated GWAS summary
statistics.  
Parameters obtained from LHC-MR, LDscores of all variants, and simulated
GWAS Z-scores are already included.  
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

In the Results/ folder can be found, for each trait, a file with
variants, p-values from PRISM, and pleiotropy annotations.

``` r
library(ggplot2)

Graph <- Example_graph(ListofTraits, ParametersTable, Trait = "B4")

print(Graph)
```

If you want to use PRISM on your own data, you can use this vignette:
\[<https://github.com/martintnr/PRISM/blob/main/vignettes/PRISM_vignette.md>\]

## Contact info

You can email me at martin.tournaire@u-paris.fr
