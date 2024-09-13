#' Graph results from the example data
#'
#' `Example_graph()` returns a ggplot of the results obtained from the example.
#'
#' @param Trait The trait of interest.
#' @param ListofTraits All traits in the analysis (necessary to calculate the p-value threshold)
#' @param ThreshSelectionPvalues P-value threshold to select PRISM top variants, if you do not want
#' the default Bonferroni correction.
#'@param ParametersTable The parameters table for each pair of traits.
#'
#' @return a ggplot of all genetic variants with GWAS p-value on the x-axis, PRISM p-value
#' on the y-axis, and colors according to pleiotropy.
#' @export
#'
#' @examples
Example_graph <- function(ListofTraits, Trait, ParametersTable, ThreshSelectionPvalues = 5e-08/(length(ListofTraits)-1)){


    Somme <- fread(paste0("Results/Pleio_", Trait,".csv"))
    path <- paste0("Data/", list.files("Data/", pattern = paste0("^",Trait,".csv")))
    X <- fread(path)

    Somme$PvalPRISM[Somme$PvalPRISM < 1e-200] <- 1e-200 #To avoid disproportionate graphs


    nX <- c(ParametersTable$nX[ParametersTable$X == Trait],
      ParametersTable$nY[ParametersTable$Y == Trait])[1]

    Somme$PvalGWAS <-  2*pnorm(q=abs(X$Zscore), lower.tail=FALSE)
    message("Pvalues from GWAS calculated from Z-scores")
    Somme$PvalGWAS[Somme$PvalGWAS < 1e-200] <- 1e-200

    Somme$SynthPleio[ Somme$SynthPleio == "No supplementary info"] <- "Direct Effect"
    Somme$SynthPleio[ Somme$SynthPleio == "Detected Network Pleiotropy"] <- "Confounder Pleiotropy"
    Somme$SynthPleio[ Somme$SynthPleio == "Suspected Vertical Pleiotropy"] <- "Vertical Pleiotropy"
    Somme$SynthPleio[ Somme$SynthPleio == "Direct Effect" & Somme$PvalPRISM > ThreshSelectionPvalues ] <- "No Effect"


    Somme$SynthPleio <-   factor( Somme$SynthPleio, levels= c("Direct Effect",
                                                             "Vertical Pleiotropy",
                                                             "Confounder Pleiotropy",
                                                             "No Effect") )


    p <- qplot(-log10(Somme$PvalGWAS), -log10(Somme$PvalPRISM), data = Somme, colour = Somme$SynthPleio,
          main = paste0("Significance of the variant-trait effects on ", Trait, ", ", nrow(Somme), " variants")
          ,  xlab = expression(paste("-",log[10] ,"(GWAS association p-value)")), ylab = expression(paste("-",log[10] ,"(PRISM p-value)"))) +
      scale_color_manual(values = c("Confounder Pleiotropy" = "#a23c33",
                                    "Direct Effect"="#45709d",
                                    "Vertical Pleiotropy"="#82992a",
                                    "No Effect" = "#D3D3D3")) +
      labs(color = "") + theme_bw() +
      theme(text = element_text(size=20))

    return(p)
  }

