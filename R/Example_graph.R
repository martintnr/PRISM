#' Graph results from the example data
#'
#' `Example_graph()` returns a ggplot of the results obtained from the example.
#'
#' @param Trait The trait of interest.
#' @param ListofTraits All traits in the analysis (necessary to calculate the p-value threshold)
#' @param ThreshSelectionPvalues P-value threshold to select PleioVar top variants, if you do not want
#' the default Bonferroni correction.
#'@param ParametersTable The parameters table for each pair of traits.
#'
#' @return a ggplot of all genetic variants with GWAS p-value on the x-axis, PleioVar p-value
#' on the y-axis, and colors according to pleiotropy.
#' @export
#'
#' @examples
Example_graph <- function(ListofTraits, Trait, ParametersTable, ThreshSelectionPvalues = 5e-08/length(ListofTraits)){


  Somme <- fread(paste0("Results/Pleio_", Trait,".csv"))
  path <- paste0("Data/", list.files("Data/", pattern = paste0("^",Trait,".csv")))
  X <- fread(path)

  Somme$PvalPleioVar[Somme$PvalPleioVar < 1e-200] <- 1e-200 #Pour éviter des graphes bizarres


  nX <- c(ParametersTable$nX[ParametersTable$X == Trait],
    ParametersTable$nY[ParametersTable$Y == Trait])[1]

  Somme$PvalGWAS <- 2*pnorm(q=abs(X$Zscore * sqrt(nX)), lower.tail=FALSE)
  message("Pvalues from GWAS calculated from Zscores and sample size")
  Somme$PvalGWAS[Somme$PvalGWAS < 1e-200] <- 1e-200

  #Faisons un graphe
  Somme$SynthPleio[ Somme$SynthPleio == "No supplementary info"] <- "\n Direct Effect \n "
  Somme$SynthPleio[ Somme$SynthPleio == "Detected Network Pleiotropy"] <- "\n Detected Network \n Pleiotropy \n "
  Somme$SynthPleio[ Somme$SynthPleio == "Suspected Vertical Pleiotropy"] <- "\n Suspected Vertical \n Pleiotropy \n "
  Somme$SynthPleio[ Somme$SynthPleio == "\n Direct Effect \n " & Somme$PvalPleioVar > ThreshSelectionPvalues ] <- "\n No Effect \n "



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

