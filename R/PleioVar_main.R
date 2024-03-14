#' Title
#'
#' @param ListofTraits
#' @param ParametersTable
#' @param Index
#' @param gzip
#' @param pU
#' @param TreshSelectionPvalues
#' @param NbCores
#' @param keepIntermediateFiles
#'
#' @return
#' @export
#'
#' @examples
#'
#'
PleioVar_main <- function(ListofTraits, ParametersTable, Index, NbCores = 1, gzip = F,
                          pU = 1e-05, TreshSelectionPvalues = 5e-08/length(ListofTraits),
                          keepIntermediateFiles = F){


  if(!file.exists("Pairwise/")){system("mkdir Pairwise")}
  if(!file.exists("Traitwise/")){system("mkdir Traitwise")}
  if(!file.exists("Results/")){system("mkdir Results")}


  if(length(list.files("Pairwise/")) != 0 | length(list.files("Traitwise/")) != 0){
    print("Error: the Pairwise/ and Traitwise/ folders must be empty")
    return()
  }
  # We set up the folder architecture


  # PleioVar is separated in two parts
  message("Executing pairwise pipeline...")


  Pairwise_pipeline(ListofTraits, ParametersTable, Index, NbCores, gzip, pU)

  message("Pairwise pipeline was successful")

  message("Executing traitwise pipeline...")

  Traitwise_pipeline(ListofTraits, ParametersTable, Index, NbCores, gzip, pU, TreshSelectionPvalues)

  message("Traitwise pipeline was successful")

  if(keepIntermediateFiles == F){
    unlink("Pairwise/", recursive = TRUE) # will delete directory
    unlink("Traitwise/", recursive = TRUE)

    message("Intermediate files removed")
  }




}
