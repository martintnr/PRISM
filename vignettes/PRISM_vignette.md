PRISM vignette
================

### Before we start

This vignette is a tutorial to use PRISM on your own data.  
Before starting, a few points:

- PRISM relies on LHC-MR results as input parameters. You will have to
  process your pairs of traits through LHC-MR before starting PRISM.  
- LHC-MR is limited to HapMap3 variants, and so is PRISM.
- LHC-MR and PRISM are computationally intensive methods, and
  exponentially so the more traits you add. Parallel computing is highly
  recommended, but requires substantial RAM resources.

### Installation

First, let’s install devtools and PRISM (if necessary).

``` r
if(!"devtools" %in% installed.packages()[, "Package"]) {
  install.packages("devtools") }

  if(!"PRISM" %in% installed.packages()[, "Package"]) {
  devtools::install_github("martintnr/PRISM") }
```

Then, let’s install LHC-MR from <https://github.com/LizaDarrous/lhcMR>
and its dependencies.

``` r
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("GenomicSEM/GenomicSEM")

devtools::install_github("LizaDarrous/lhcMR")
```

### Setting up the data

Let’s create and use as working directory PRISM_NewData folder that will
be used to store necessary and temporary files.

``` r

if(!grepl("PRISM_New", getwd(), fixed = TRUE)){
if(!file.exists("PRISM_New/")){system("mkdir PRISM_New")}
setwd("PRISM_New")}
```

Some data (like LD-scores) are required for PRISM and its dependencies.
Let’s download them from zenodo, and also create a folder for the future
LHC-MR results.

``` r

system("wget -O Necessary_data.zip  https://zenodo.org/records/10829652/files/Necessary_data.zip?download=1")
system("unzip Necessary_data.zip -d Necessary_data && rm Necessary_data.zip")


if(!file.exists("LHCMR_Results")){system("mkdir LHCMR_Results")}
```

Next, you have to indicate the path of your GWAS data and the traits you
want to process.  
For this example, I will use LDL, HDL, and CAD from UKBiobank round 2.  
You should use LHC-MR (<https://github.com/LizaDarrous/lhcMR>)
guidelines to format your data.  
Please format your GWAS data filenames to begin with *Trait*. (the
default UKBiobank round 2 format works). If multiple files in the folder
begins with the same prefix, please remove or rename the duplicates.  
Both .tsv and .csv formats are acceptable, and the data can be gzipped
(.gz or .bgz).

``` r

#system("mkdir GWAS_data")

#system("wget -O HDL.tsv.bgz  https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30760_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz")


#system("wget -O LDL.tsv.bgz  https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30780_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz")

#system("wget -O CAD.tsv.bgz  https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/I25.gwas.imputed_v3.both_sexes.tsv.bgz")

#system("mv *.bgz  GWAS_data/")
```

### Running LHC-MR on all pairs of traits

First, please put the path to the folder containing your data in the
`sourceGWAS` object, and a `ListeofTraits` vector with all your
traits.  
The `AllPairs` object will contain all possible unique pairs of traits.

``` r
# sourceGWAS <- "GWAS_data/"
# ListofTraits <- c("LDL", "HDL", "CAD")
AllPairs <- (t(combn(ListofTraits,2)))
colnames(AllPairs) <- c("X", "Y")
```

Now, you can process your traits with LHC-MR. The function below is
simply a loop of the example from their github.  
The `NbCores` parameter is passed to `lhc_mr()` and is very important to
speed up calculation. Beware, a high number of cores will greatly speed
up the process, but will also consume a lot of RAM.  
Results will be written in the `LHCMR_Results`folder.

``` r
  
  library(data.table)
  library(GenomicSEM)
  library(TwoSampleMR)
  library(lhcMR)

  RefFolder <- getwd()

  VAR = fread(paste0(RefFolder, "/Necessary_data/variants.tsv.bgz"),
            select = c("variant", "chr", "pos", "ref", "alt", "rsid"))
  VAR$chr=as.numeric(VAR$chr)

lhcmrLoop <- function(A, NbCores, Minimum_MAF = 0.05, run_ldsc = T, run_MR = T){
  gc()
  print(paste0("Loop"," ",A)) 
  
    setwd(RefFolder)

    debut <- Sys.time()

  Trait1 <- AllPairs[A,"X"] # Get the trait from the specific pairing
  Trait2 <- AllPairs[A,"Y"]

  

  # Import both sets of summary stats and remove missing observations (= variants with no calculable P-values)
  
  path1 <- paste0(sourceGWAS, list.files(sourceGWAS, pattern = paste0("^",Trait1)))
  X <- fread(path1)
 
  path2 <- paste0(sourceGWAS, list.files(sourceGWAS, pattern = paste0("^",Trait2)))
  Y <- fread(path2)
 


  
  # Low confidence variants and variants with low minor allele frequencies (MAF) are filtered out.
  X <- X[ (! is.na(X$pval)) & (! X$low_confidence_variant) & ( X$minor_AF > Minimum_MAF), ]
  Y <- Y[ (! is.na(Y$pval)) & (! Y$low_confidence_variant) & ( Y$minor_AF > Minimum_MAF), ] 
  #
  
  
  ## To make sure we have all the columns needed (including allele data), we merge with Neale-provided variants file

  X = dplyr::inner_join(X,VAR[,c(1:6)])
  Y = dplyr::inner_join(Y,VAR[,c(1:6)])
  
  ## File paths needed for the analysis
  LD.filepath = paste0(RefFolder, "/Necessary_data/LDscores_filtered.csv") # LD scores
  rho.filepath = paste0(RefFolder, "/Necessary_data/LD_GM2_2prm.csv") # local/SNP-specific LD scores
  
  
  
  
  
  ld = paste0(RefFolder, "/Necessary_data/eur_w_ld_chr/")  #LD information
  hm3 = paste0(RefFolder, "/Necessary_data/w_hm3.snplist")
  
 
  ## Step 1
  trait.names=c(Trait1,Trait2)
  input.files = list(X,Y)
  
  setwd("LHCMR_Results") #LHC-MR saves a lot of stuff in the working directory
  
  df = merge_sumstats(input.files,trait.names,LD.filepath,rho.filepath) #code from LHCMR
  #save(df, file = "gen_data")
  rm(X)
  rm(Y)
  rm(input.files)
  gc()
  
  ## Step 2
  SP_list = calculate_SP(df,trait.names,run_ldsc,run_MR,hm3=hm3,ld=ld,nStep = 2,
                         SP_single=3,SP_pair=100,SNP_filter=10, nCores = NbCores) #Calculating the starting points
  
  ## Step 3
  res = lhc_mr(SP_list, trait.names, paral_method= "lapply", nBlock=200, 
               nCores = NbCores, run_ldsc, run_MR) #The actual optimisation
  
  res <- as.data.frame(res)
  
  res$nX <- c(unique(df$`N.x`), NA, NA)
  res$nY <- c(unique(df$`N.y`), NA, NA)

  Nom = paste0("ResultsLHCMR_", Trait1, "_", Trait2)
  save(res, file = Nom)
  
  
  fin <- Sys.time()
  print(c(Trait1,Trait2))
  print(fin-debut)
  
  setwd(RefFolder)

  
  return(res)
}


Range <- c(1:nrow(AllPairs))  

lapply(X = Range, FUN = lhcmrLoop, NbCores = 1, run_MR = F)
```

### Running PRISM

Let’s format the LHC-MR results in a PRISM-friendly format, and gather
the last necessary inputs.

``` r

library(data.table)
library(dplyr)
library(stats)
library(stringr)
library(parallel)
library(matrixStats)
library(PRISM)



ParametersTable <- Format_parameters_table(AllPairs, resdir = "LHCMR_Results")
Index <- fread("Necessary_data/Index_FULL.csv", header = T, sep = ",")
colnames(Index) <- c("rsid", "variant", "LDscore")
ListofTraits <- unique(c(ParametersTable$X, ParametersTable$Y))

  
```

The final step is to run the PRISM pipeline.  
If you have a *standard_beta* column in your GWAS data, PRISM will
process them.  
If you do not, the standardized effect sizes will be computed with
$\frac{tstat}{\sqrt{(sample \ size)}}$. If you formatted your data
according to LHC-MR guidelines and the previous steps worked, everything
should run automatically.  
As usual, a higher number of cores in the `NbCores` parameters will
speed up calculation.

``` r

PRISM_main(ListofTraits, ParametersTable, Index , sourceGWAS = sourceGWAS,
              NbCores = 1, gzip = F, keepIntermediateFiles = F)
```

### [Pre-print](https://doi.org/10.1101/2024.06.01.24308193) results

The results presented in our pre-print can be reproduced using LHC-MR
results available [here](https://doi.org/10.5281/zenodo.13862821) and UK
Biobank GWAS round 2 summary statistics available
[here](https://www.nealelab.is/uk-biobank). The full list of processed
traits is provided in [Supplementary Table
2](https://www.medrxiv.org/content/10.1101/2024.06.01.24308193v2.supplementary-material).
