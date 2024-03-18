#' Format parameters table
#'
#' @param resdir
#' @param AllPairs
#' @param rhoXY
#'
#' @return
#' @export
#'
#' @examples
Format_parameters_table <- function(AllPairs, resdir, rhoXY = 0){


  colnames(AllPairs) <- c("X", "Y")

  get_parameters <- function(A){

    Trait1 <- AllPairs[A,"X"]
    Trait2 <- AllPairs[A,"Y"]

    Nom = paste0(resdir,"/ResultsLHCMR_", Trait1, "_", Trait2)

    load(Nom) #res object is loaded

    params <- c(Trait1, Trait2,
                res[1,"iX"], res[1,"piX"],res[1,"h2X"],res[1,"tX"],res[1,"axy"],res[3,"axy"],res[1,"nX"],
                res[1,"iY"], res[1,"piY"],res[1,"h2Y"],res[1,"tY"],res[1,"ayx"],res[3,"ayx"],res[1,"nY"],
                rhoXY)

    return(params)
  }
  tabparams <- as.data.frame(do.call(rbind, lapply(X = c(1:nrow(AllPairs)), FUN =  get_parameters)))


  colnames(tabparams) <- c("X", "Y",
                           "iX", "piX", "h2X", "tX", "axy", "pval_axy","nX",
                           "iY", "piY", "h2Y", "tY", "ayx", "pval_ayx", "nY",
                           "rhoXY")

  return(tabparams)

}


