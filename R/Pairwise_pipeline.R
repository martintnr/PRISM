#' Title
#'
#' @param ListofTraits
#' @param ParametersTable
#' @param NbCores
#' @param Index
#' @param gzip
#' @param pU
#'
#' @return
#' @export
#'
#' @examples
#'
#'
Pairwise_pipeline <- function(ListofTraits, ParametersTable, Index ,NbCores, gzip, pU){




    CharcClassif <- function(A){




      gc()
      Trait1 <-  ParametersTable$X[A]
      Trait2 <-  ParametersTable$Y[A]

      print(paste0("Pair ", Trait1, " ", Trait2))

      path1 <- paste0("Data/", list.files("Data/", pattern = paste0("^",Trait1,".csv")))
      path2 <- paste0("Data/", list.files("Data/", pattern = paste0("^",Trait2,".csv")))

      X <- fread(path1)
      Y <- fread(path2)



      nX <- ParametersTable$nX[A]
      nY <- ParametersTable$nY[A]
      rho <- ParametersTable$rhoXY[A]

      M = nrow(Index) #Number of SNPs
      bX = X$Zscore  #Effects of trait X
      bY = Y$Zscore  #Effects of trait Y

      ld = Index$LDscore



      ### Calculates log likelihood from parameters



      ld2 = ld^2

      iX = ParametersTable$iX[A]
      iY = ParametersTable$iY[A]
      pX = ParametersTable$piX[A]
      pY = ParametersTable$piY[A]
      h2X = ParametersTable$h2X[A]
      h2Y = ParametersTable$h2Y[A]
      tX = ParametersTable$tX[A]
      tY = ParametersTable$tY[A]
      a = ParametersTable$axy[A]
      b = ParametersTable$ayx[A]


      L = NA

      T0 = matrix(c(iX/nX, rho, iY/nY), ncol = 3, nrow = 1)
      T1 = (h2X / (M*pX)) * matrix(c(1, a, a^2), ncol = 3, nrow = 1)
      T2 = (1 / (M*pU)) * matrix(c((tX + (b*tY))^2, (tX + (b*tY)) * (tY + (a*tX)), (tY + (a*tX))^2), ncol = 3, nrow = 1)
      T3 = (h2Y / (M*pY)) * matrix(c(b^2, b, 1), ncol = 3, nrow = 1)

      detS0 = (T0[1]*T0[3])-(T0[2]^2)


      U = T1
      detA = (U[1]*U[3])-(U[2]^2)
      detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
      det0 = (ld2*detA)+(ld*detB)+detS0
      dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
      l01 = pX*(1-pU)*(1-pY)*dmv_sub #Om1

      U = T2;
      detA = (U[1]*U[3])-(U[2]^2)
      detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
      det0 = (ld2*detA)+(ld*detB)+detS0
      dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
      l02 = pU*(1-pX)*(1-pY)*dmv_sub #Om2

      U = T3;
      detA = (U[1]*U[3])-(U[2]^2)
      detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
      det0 = (ld2*detA)+(ld*detB)+detS0
      dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
      l03 = pY*(1-pX)*(1-pU)*dmv_sub #om3

      U = T1+T2;
      detA = (U[1]*U[3])-(U[2]^2)
      detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
      det0 = (ld2*detA)+(ld*detB)+detS0
      dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
      l012 = pX*pU*(1-pY)*dmv_sub #Om4

      U = T1+T3;
      detA = (U[1]*U[3])-(U[2]^2)
      detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
      det0 = (ld2*detA)+(ld*detB)+detS0
      dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
      l013 = pX*pY*(1-pU)*dmv_sub #Om5

      U = T2+T3;
      detA = (U[1]*U[3])-(U[2]^2)
      detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
      det0 = (ld2*detA)+(ld*detB)+detS0
      dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
      l023 = pU*pY*(1-pX)*dmv_sub #Om6

      U = T1+T2+T3;
      detA = (U[1]*U[3])-(U[2]^2)
      detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
      det0 = (ld2*detA)+(ld*detB)+detS0
      dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
      l0123 = pX*pU*pY*dmv_sub #Om7

      l0 = (1-pX)*(1-pU)*(1-pY)*(detS0^(-1/2))*(exp( (-1/2)*(1/detS0) * ((T0[3]*(bX^2)) + (T0[1]*(bY^2)) - (2*T0[2]*(bX*bY))))); #Om0

      L = log(l01+l02+l03+l012+l013+l023+l0123+l0)-log(2*pi)
      L = -sum(L) #unpenalized global likelihood, not really used



      #Likelihood of each variant to belong to each model
      Om1 = l01 # 3+2+2 = 7 parameters
      Om2 = l02 # 3+4+2 = 9 parameters
      Om3 = l03 # 3+2+2 = 7 parameters
      Om4 = l012 #3+2+4+2 = 11
      Om5 = l013 #3+2+2+2 = 9
      Om6 = l023 #3+2+4+2 = 11
      Om7 = l0123  #3+2+4+2+2 = 13
      Om0 = l0  #3+2 = 5 parameters

      Omega.list <- list(L,Om0,Om1,Om2,Om3,Om4,Om5,Om6,Om7)

      vraisembl <- unlist(lapply(Omega.list, '[[', 1)) #On récupère les vraisemblances
      indice <- which.min(vraisembl)
      Omega_opti <- as.data.frame(Omega.list)
      colnames(Omega_opti) = c("L2", "Om0", "Om1", "Om2", "Om3", "Om4", "Om5", "Om6", "Om7") #On a récupérer la vraisemblance pour chaque modèle pour chaque SNP


      Omega_opti$INDEX <- Index$variant


      write.table(Omega_opti, paste0("Pairwise/Likelihood_", Trait1, "_", Trait2,".csv"), sep=",", quote=F, row.names=F)

      if(gzip == T){
        tryCatch(
          { system(paste0("gzip Pairwise/Likelihood_", Trait1, "_", Trait2, ".csv"))


          },
          error = function(cond) {
            message("gzip did not work but that's okay")
          })
      }


      M = nrow(Omega_opti)


      Omega_BIC <- Omega_opti
      Omega_BIC <- transform(Omega_BIC, index =1+apply(Omega_BIC[,2:9], 1, which.max))
      Label <- c( "L2", "Om0", "Om1", "Om2", "Om3", "Om4", "Om5", "Om6", "Om7")
      if(!"SNP_Category"%in% colnames(Omega_BIC)){ #On ne remplit le tableau avec la catégorie que si cela n'a pas été fait dans un des if précédent
        Omega_BIC$SNP_Category <- Label[Omega_BIC$index]}else if(is.na(Omega_BIC$SNP_Category[1])){
          Omega_BIC$SNP_Category <- Label[Omega_BIC$index]
        }
      freqs <-t(apply(Omega_BIC[,2:9] ,1, function(x) x/sum(x)))
      Score = data.frame(PH = 1:nrow(Omega_BIC)) #PlaceHolder
      Score[[paste0("O_",Trait1)]] <- freqs[,"Om0"] #Effet sur rien
      Score[[paste0("O_",Trait2)]] <- freqs[,"Om0"] #Effet sur rien
      Score[[paste0("E_",Trait1)]] <- apply(freqs[,c(2,5,6,8)], 1, max, na.rm=TRUE) #Effet sur X
      Score[[paste0("U_",Trait1,"_",Trait2)]] <- apply(freqs[,c(3,5,7,8)], 1, max, na.rm=TRUE) #Effet sur U
      Score[[paste0("E_",Trait2)]] <-apply(freqs[,c(4,6,7,8)], 1, max, na.rm=TRUE) #Effet sur Y
      Score$PH <- NULL
      Score$variant <- Omega_BIC$INDEX
      Score$CatCharc <- Omega_BIC$SNP_Category


      write.table(Score, paste0("Pairwise/SNP_Scores_", Trait1, "_", Trait2,".csv"), sep=",", quote=F, row.names=F, col.names = T)


      if(gzip == T){
        tryCatch(
          { system(paste0("gzip Pairwise/SNP_Scores_", Trait1, "_", Trait2,".csv"))


          },
          error = function(cond) {
            message("gzip did not work but that's okay")
          })
      }


    }

    RANGE <- c(1:nrow(ParametersTable))
    mclapply(X = RANGE, FUN = CharcClassif, mc.cores = NbCores)

}
