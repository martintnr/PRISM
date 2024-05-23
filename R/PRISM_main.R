#' Main pipeline of PRISM
#' @description
#' `PRISM_main()` is the main function of PRISM, that takes as input the list of traits to process,
#' the LDscore index, and the LHC-MR parameters table.
#'
#' It will write in the Results/ folder, for each trait, the pleiotropic label
#' of all variants.
#'
#'
#' @param ListofTraits The list of traits to be processed by PRISM
#' @param ParametersTable The parameters table for each pair of traits, obtained
#' from LHC-MR
#' @param Index A dataframe with 2 columns, that contain respectively the name of
#' the variants and their LDscores.
#' @param gzip TRUE if you want the example data to be gzipped, FALSE otherwise.
#' gzip= F will speed up computations, but will increase disk space usage (only
#' relevant if keepIntermediateFiles = TRUE).
#' @param pU Polygenicity of the confounder U, as this parameter is not obtained
#' from LHC-MR. 1e-05 is a robust value, we do not advise to change it.
#' @param ThreshSelectionPvalues P-value threshold to select PRISM top variants, if you do not want
#' the default Bonferroni correction.
#' @param NbCores The number of cores to use, passed to `mclapply()`. Greatly
#' reduces the computation time.
#' @param keepIntermediateFiles TRUE if you want to keep the intermediates files,
#' which can take a lot of disk space with a high number of traits and genetic variants.
#' We advise to set gzip = T in this case. FALSE if you want to remove them.
#' @param sourceGWAS The folder where GWAS summary statistics files are held.
#' @param Minimum_MAF Genetic variants with minimum allele frequency below this value will be filtered out.
#' @param skipPairwise TRUE if you only want to execute the traitwise pipeline.
#' @param labelGWASsig TRUE if you want to also label variants that are not PRISM significant but are GWAS significant with the same threshold.
#' A "variant" and a "pval" columns are necessary in the GWAS summary statistics files.
#'
#' @return
#' @export
#'
#' @examples
#'
#'
PRISM_main <- function(ListofTraits, ParametersTable, Index, sourceGWAS = getwd() ,NbCores = 1, gzip = F,
                          pU = 1e-05, ThreshSelectionPvalues = 5e-08/length(ListofTraits),
                          keepIntermediateFiles = F,  Minimum_MAF = 0.05, skipPairwise = F, labelGWASsig = F){


  if(!file.exists("Pairwise/")){system("mkdir Pairwise")}
  if(!file.exists("Traitwise/")){system("mkdir Traitwise")}
  if(!file.exists("Results/")){system("mkdir Results")}


  if(length(list.files("Pairwise/")) != 0 | length(list.files("Traitwise/")) != 0){
    print("Error: the Pairwise/ and Traitwise/ folders should be empty")
  #  return()
  }
  # We set up the folder architecture


  if(length(ListofTraits) < 31){
    message(paste0(length(ListofTraits)), " traits were submitted, which is less than the 31 traits recommended. Power and accuracy will be reduced")
    }

  # PRISM is separated in two parts

  if(skipPairwise == F){
  message("Executing pairwise pipeline...")


  Pairwise_pipeline(ListofTraits, ParametersTable, Index , sourceGWAS, NbCores, gzip, pU,  Minimum_MAF)

  message("Pairwise pipeline was successful")

  }else{  message("Pairwise pipeline was skipped, as instructed")}

  message("Executing traitwise pipeline...")

  Traitwise_pipeline(ListofTraits, ParametersTable, Index ,sourceGWAS, NbCores, gzip, pU, ThreshSelectionPvalues, labelGWASsig)

  message("Traitwise pipeline was successful")

  if(keepIntermediateFiles == F){
    unlink("Pairwise/", recursive = TRUE) # will delete directory
    unlink("Traitwise/", recursive = TRUE)

    message("Intermediate files removed")
  }




}
