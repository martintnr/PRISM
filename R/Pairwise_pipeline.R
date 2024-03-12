#' Title
#'
#' @param ListofTraits
#' @param ParametersTable
#' @param NbCores
#'
#' @return
#' @export
#'
#' @examples
#'
#'
Pairwise_pipeline <- function(ListofTraits, ParametersTable, NbCores){


  Combinaison <- (t(combn(ListofTraits,2))) #All the pairs


    CharcClassif <- function(A){




      gc()
      Trait1 <-  Combinaison[A,1]
      Trait2 <-  Combinaison[A,2]

      exp <- fread()
      out <- fread()





    }

}
