#' Title
#'
#' @param Trait
#' @param ListofTraits
#' @param TreshSelectionPvalues
#'
#' @return
#' @export
#'
#' @examples
Example_graph <- function(ListofTraits, Trait, TreshSelectionPvalues = 5e-08/length(ListofTraits)){


  Somme <- fread(paste0("Results/Pleio_", Trait,".csv"))
  path <- paste0("Data/", list.files("Data/", pattern = paste0("^",Trait,".csv")))
  X <- fread(path)

  Somme$PvalPleioVar[Somme$PvalPleioVar < 1e-200] <- 1e-200 #Pour éviter des graphes bizarres

  Somme$PvalGWAS <- 2*pnorm(q=abs(X$Zscore/0.005), lower.tail=FALSE)
  message("Pvalues from GWAS calculated from from ZScores")
  Somme$PvalGWAS[Somme$PvalGWAS < 1e-200] <- 1e-200

  #Faisons un graphe
  Somme$SynthPleio[ Somme$SynthPleio == "No supplementary info"] <- "\n Direct Effect \n "
  Somme$SynthPleio[ Somme$SynthPleio == "Detected Network Pleiotropy"] <- "\n Detected Network \n Pleiotropy \n "
  Somme$SynthPleio[ Somme$SynthPleio == "Suspected Vertical Pleiotropy"] <- "\n Suspected Vertical \n Pleiotropy \n "
  Somme$SynthPleio[ Somme$SynthPleio == "\n Direct Effect \n " & Somme$PvalPleioVar > TreshSelectionPvalues ] <- "\n No Effect \n "



  p <- qplot(-log10(Somme$PvalGWAS), -log10(Somme$PvalPleioVar), data = Somme, colour = Somme$SynthPleio,
        main = paste0("Significance of the effect of variants on ", Trait, ", ", nrow(Somme), " variants")
        , xlab = "-log10(GWAS association p-value)", ylab = "-log10(PleioVar p-value)") +
    scale_color_manual(values = c("\n Detected Network \n Pleiotropy \n " = "#a23c33",
                                  "\n Direct Effect \n "="#45709d",
                                  "\n Suspected Vertical \n Pleiotropy \n "="#82992a",
                                  "\n No Effect \n " = "#D3D3D3")) +
    labs(color = "Pleiotropy") +
    theme(text = element_text(size=20))

  return(p)
  }

