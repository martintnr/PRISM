library(data.table)

setwd("/home/PleioMap/Data/LHCMR/PleioVar_Scenarios/")


Traits <-rev(c("A", #Anciens traits triés
               "B1",
               "B2",
               "B3",
               "B4", #Our parangon
               "C1",
               "C2",
               "C3",
               "C4",
               "D",
               "E1",
               "E2",
               "E3",
               "E4",
               "E5",
               "E0"))


Combinaison <- t(combn(Traits,2))


Scenario_parameters <- read.csv(file = "GWAS_Simulated/Scenario_parameters.csv", header = T)
Causal_Relationships <- read.csv(file = "GWAS_Simulated/Causal_Relationships.csv", header = T)


### Get all simulated data
GetAllData <- function(A, ID = 56, replica = 1){

  Trait <- A
  exp <- fread(paste0("GWAS_Simulated/SimuDataFull_", ID,"_", Trait,"_", replica, ".csv.gz"))
  return(exp[,c(1:5)])
}
listedata <- lapply(X = Traits, FUN = GetAllData)



GetAllParams <- function(A, ID = 56, replica = 1){


  Trait1 <-  Combinaison[A,1]
  Trait2 <-  Combinaison[A,2]



 # exp <- fread(paste0("GWAS_Simulated/SimuDataFull_", ID,"_", Trait1,"_", replica, ".csv.gz"))


  #out <- fread(paste0("GWAS_Simulated/SimuDataFull_", ID,"_", Trait2,"_", replica, ".csv.gz"))




  nX=361194	# Sample size of the exposure (celle de UK biobank )
  nY=361194





  ### Read in data from previously generated data
  #M = nrow(exp) #Number of SNPs
  phen_corr = 0
  nXY = (nX/380e3)*(nY/380e3)*380e3 #Estimated fraction of sample overlap
  rho = phen_corr*nXY/(as.numeric(nX)*as.numeric(nY)) #Calculated phenotypic correlation due to sample overlap, ici vaut 0


  iX = 1.2
  iY = 1.2
  pX = Scenario_parameters$Polygenicity[ID]
  pY = Scenario_parameters$Polygenicity[ID]

  if(Trait1 == "A" | Trait1 == "B1" |Trait1 == "B2" |Trait1 == "B3" |Trait1 == "B4" ){h2X <- Scenario_parameters$h2_AB[ID]
  }else{h2X <- Scenario_parameters$h2_CDE[ID]}

  if(Trait2 == "A" | Trait2 == "B1" |Trait2 == "B2" |Trait2 == "B3" |Trait2 == "B4" ){h2Y <- Scenario_parameters$h2_AB[ID]
  }else{h2Y <- Scenario_parameters$h2_CDE[ID]}

  Check = Causal_Relationships[Causal_Relationships$Traits == Trait2, Trait1]
  if(Check == "Causal effect" & Trait1 == "B4"){b <- Scenario_parameters$Causal_toB4[ID]}
  if(Check == "Causal effect" & Trait1 != "B4"){b <- Scenario_parameters$Causal_exceptB4[ID]}
  if(Check == "No causal effect"){b <- 0}

  Check = Causal_Relationships[Causal_Relationships$Traits == Trait1, Trait2]
  if(Check == "Causal effect" & Trait2 == "B4"){a <- Scenario_parameters$Causal_toB4[ID]}
  if(Check == "Causal effect" & Trait2 != "B4"){a <- Scenario_parameters$Causal_exceptB4[ID]}
  if(Check == "No causal effect"){a <- 0}

  tX = 0.05
  tY = 0.05


  pU  <-  Scenario_parameters$pU_charcut[ID]


  pvala <- 1e-10
  pvalb <- 1e-10

  if(a == 0){pvala <- 1}
  if(b == 0){pvalb <- 1}




  return(  c(list(Trait1, Trait2, iX, pX, h2X, tX, a, pvala, iY, pY, h2Y, tY, b, pvalb)))

}

listeparams <- lapply(X = c(1:nrow(Combinaison)), FUN = GetAllParams)


tabparams <- as.data.frame(do.call(rbind, listeparams))
colnames(tabparams) <- c("X", "Y", "iX", "piX", "h2X", "tX", "axy", "pval_axy",
                         "iY", "piY", "h2Y", "tY", "ayx", "pval_ayx")

tabparams$X <- unlist(tabparams$X)
tabparams$Y <- unlist(tabparams$Y)

tabparams[,3:ncol(tabparams)] <-  as.numeric(unlist(tabparams[,3:ncol(tabparams)]))

simulated_example_data <- list(listedata, tabparams)
usethis::use_data(simulated_example_data)


##### Et l'index



VAR2 = fread( "/home/PleioMap/Data/LHCMR/PleioVarNext2/Likelihood/VAR_LIGHT")

usethis::use_data(VAR2)

