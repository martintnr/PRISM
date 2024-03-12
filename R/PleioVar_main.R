#' Title
#'
#' @param ListofTraits
#' @param ParametersTable
#' @param NbCores
#' @return
#' @export
#'
#' @examples
#'
#'
PleioVar_main <- function(ListofTraits, ParametersTable, NbCores){


  system("mkdir Pairwise")
  system("mkdir Traitwise")


Pairwise_pipeline(ListofTraits, ParametersTable, NbCores)
Traitwise_pipeline(ListofTraits, ParametersTable, NbCores)


}
