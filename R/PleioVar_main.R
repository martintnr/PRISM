#' Title
#'
#' @param ListofTraits
#' @param ParametersTable
#' @param Index
#' @param gzip
#' @param pU
#' @param TreshSelectionPvalues
#' @param NbCores
#'
#' @return
#' @export
#'
#' @examples
#'
#'
PleioVar_main <- function(ListofTraits, ParametersTable, Index, NbCores = 1, gzip = F,
                          pU = 1e-05, TreshSelectionPvalues = 5e-08/length(ListofTraits)){


  if(!file.exists("Pairwise/")){system("mkdir Pairwise")}
  if(!file.exists("Traitwise/")){system("mkdir Traitwise")}
  if(!file.exists("Results/")){system("mkdir Results")}


  if(length(list.files("Pairwise/")) != 0 | length(list.files("Traitwise/")) != 0){
    print("Error: the Pairwise/ and Traitwise/ folders must be empty")
    return()
  }
  # We set up the folder architecture


  # PleioVar is separated in two parts
  Pairwise_pipeline(ListofTraits, ParametersTable, Index, NbCores, gzip, pU)

  message("Pairwise pipeline was successful")

  Traitwise_pipeline(ListofTraits, ParametersTable, Index, NbCores, gzip, pU, TreshSelectionPvalues)

  message("Traitwise pipeline was successful")







}
