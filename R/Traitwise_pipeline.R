#' Traitwise pipeline of PRISM
#' @description
#' `Traitwise_pipeline()` handles the traitwise process of PRISM. It takes as input
#' the likelihoods and scores obtained from `Pairwise_pipeline()`. It outputs the syntheses
#' and the PRISM p-values by traits in the Traitwise/ folder.
#'
#' Final results (PRISM p-values and pleiotropic labels, for all variants) are
#' written for each trait in the Results/ folder.
#'
#' All parameters are passed from `PRISM_main()` and are more detailed over there.
#'
#' @param ListofTraits The list of traits to be processed.
#' @param ParametersTable The parameters table for each pair of traits.
#' @param NbCores The number of cores to use.
#' @param Index A dataframe with all variants and their LDscores.
#' @param gzip TRUE to gzip created files, FALSE otherwise.
#' @param pU Polygenicity of the confounder U.
#' @param ThreshSelectionPvalues P-value threshold to select PRISM top variants.
#' @param sourceGWAS The folder where GWAS summary statistics files are held.
#' @param labelGWASsig TRUE if you want to label variants that are not PRISM significant but are GWAS significant.
#'
#' @return
#' @export
#'
#' @examples
Traitwise_pipeline <- function(ListofTraits, ParametersTable, Index ,sourceGWAS, NbCores, gzip, pU, ThreshSelectionPvalues, labelGWASsig, Test){

  TableBase <- c(Index$variant)


  SyntheseTraits <- function(X){ #Just used to regroup data by traits for post-processing


    # Synthese
    TRAIT <- X

    if(file.exists( paste0("Traitwise/Pvalues_", TRAIT, ".csv")) | file.exists(paste0("Traitwise/Pvalues_", TRAIT, ".csv.gz"))){
      print(paste0("Pvalues for", TRAIT, " already exist. Skipped."))
      return(NA)
    }

    print(paste0(TRAIT, " synthesis"))
  #  print(paste0("gzip is ", gzip))

    gc()
    TableBaseOEffect <- as_tibble(TableBase)
    TableBaseXEffect <- as_tibble(TableBase)
    TableBaseUEffect <- as_tibble(TableBase)

    TableRecapOm2 <- as_tibble(TableBase)



    file.ls <- list.files(path="Pairwise", pattern=paste0("_",TRAIT, "_"))
    file.ls <- c(file.ls, list.files(path="Pairwise", pattern=paste0("_",TRAIT, "\\.csv")) )

    file.ls <- grep(paste0("^SNP_Scores_"), file.ls, value = TRUE)

    # We get all the files of the TRAIT


    for(FileNumber in c(1:length(file.ls))){

      gc()
   #    print(FileNumber)
      Scores <- fread(paste0("Pairwise/",file.ls[FileNumber]))

      TableRecapOm2 <- cbind(TableRecapOm2, Scores$CatCharc)

      Scores <- Scores %>% select(variant,contains(TRAIT))


      TableBaseOEffect <- cbind(TableBaseOEffect, Scores[[paste0("O_", TRAIT)]])
      TableBaseXEffect <- cbind(TableBaseXEffect, Scores[[paste0("E_", TRAIT)]])
      NomColonneU <- str_subset(colnames(Scores), "U_" )
      TableBaseUEffect <- cbind(TableBaseUEffect, Scores[[NomColonneU]])
      rm(Scores)
      gc()


    }
 #   print("Synth done")

    file.ls <- gsub(pattern = paste0("(.*SNP_Scores_)(.*)(.csv.*)"), replacement = "\\2",  x = file.ls)
    file.ls <- str_remove_all(file.ls, paste0(TRAIT,"_"))
    file.ls <- str_remove_all(file.ls, paste0("_",TRAIT))

    names(TableBaseOEffect) <- c("variant", file.ls )
    names(TableBaseXEffect) <- c( "variant", file.ls )
    names(TableBaseUEffect) <- c("variant", file.ls )
    names(TableRecapOm2) <- c("variant", file.ls )




    write.table(TableBaseOEffect, file = paste0("Traitwise/Synthese_OEffect_", TRAIT,".csv"), sep=",", quote=F, row.names=F, col.names = T)
    write.table(TableBaseXEffect, file = paste0("Traitwise/Synthese_XEffect_", TRAIT,".csv"), sep=",", quote=F, row.names=F, col.names = T)
    write.table(TableBaseUEffect, file = paste0("Traitwise/Synthese_UEffect_", TRAIT,".csv"), sep=",", quote=F, row.names=F, col.names = T)
    write.table(TableRecapOm2, file = paste0("Traitwise/Synthese_CatCharc_", TRAIT,".csv"), sep=",", quote=F, row.names=F, col.names = T)

 #   print("Start gzip S")
#    print(gzip)

    if(gzip == T){
  #    print("gzip loop entered")
      tryCatch(
        {     system(paste0("gzip Traitwise/Synthese_*", "_",TRAIT,".csv"))

        },
        error = function(cond) {
          message("gzip did not work but that's okay")
        })
    }

    ### Pre-analyse

    gc()


    TableBaseOEffect <- TableBaseOEffect[complete.cases(TableBaseOEffect), ]
    TableBaseXEffect <- TableBaseXEffect[complete.cases(TableBaseXEffect), ]

    if(nrow(TableBaseOEffect ) != nrow(TableBaseXEffect)){

      TableBaseOEffect <- TableBaseOEffect[TableBaseOEffect$variant %in% TableBaseXEffect$variant,]
      TableBaseXEffect <- TableBaseXEffect[TableBaseXEffect$variant %in% TableBaseOEffect$variant,]

    }

    gc()

    Range <- c(2:as.numeric(length(TableBaseOEffect)))


    MSD <- as.data.frame(TableBaseOEffect[,c(1)])
    colnames(MSD) <- "variant"
    MSD$meanO <-  rowMeans(TableBaseOEffect[, Range])
    MSD$sdO <- rowSds(as.matrix(TableBaseOEffect[, Range]), na.rm = TRUE)
    MSD$meanX <-  rowMeans(TableBaseXEffect[, Range])
    MSD$sdX <- rowSds(as.matrix(TableBaseXEffect[, Range]), na.rm = TRUE)
#    MSD$meanU <-  rowMeans(TableBaseUEffect[, Range])
 #   MSD$sdU <- rowSds(as.matrix(TableBaseUEffect[, Range]), na.rm = TRUE)
   # print("Start pval")


    TEST <- function(SNP){ #Paired T-test of the direct effect of the variant on the trait

      if(sum(c(na.omit(as.numeric(TableBaseXEffect[SNP,  Range])) - na.omit(as.numeric(TableBaseOEffect[SNP, Range]))) > 0.999999999999) == length(Range)){
        P = 1e-300
        #print(SNP)
        return(P)
      }


      student <-tryCatch( {t.test(as.numeric(TableBaseXEffect[SNP,  Range]), as.numeric(TableBaseOEffect[SNP, Range]),
                                  alternative= "greater", paired = TRUE)},

                          error=function(cond) {#print(SNP)

                            return( wilcox.test(as.numeric(TableBaseXEffect[SNP,  Range]), as.numeric(TableBaseOEffect[SNP, Range]),
                                                alternative= "greater", paired = TRUE))
                          }) #On compare les moyennes de nos traits, student normalement, wilcox si student crash

      P <- student$p.value
      return(P)
    }

    SIGNTEST <- function(SNP){ #Paired sign test of the effect of the variant on trait X

      leSNPX <- as.numeric(TableBaseXEffect[SNP,  Range])
      leSNPO <- as.numeric(TableBaseOEffect[SNP,  Range])
      positive_count <- sum(leSNPX > leSNPO)
      total <- positive_count + sum(leSNPX <= leSNPO)
      return(binom.test(positive_count, total, p = 0.5, alternative = "greater")$p.value)
    }


    if(Test == "Student"){ResX <- mclapply(X = c(1:nrow(MSD)), FUN = TEST, mc.cores = 1)}
    if(Test == "Sign"){ResX <- mclapply(X = c(1:nrow(MSD)), FUN = SIGNTEST, mc.cores = 1)}



    MSD$PX <- as.numeric(ResX)
    #print("end pval")


    MSD_full <-  merge(Index, MSD, by.x = "variant", by.y = "variant", no.dups = TRUE, all.x = TRUE, sort = FALSE )
    MSD_full <-  MSD_full[match(Index$variant, MSD_full$variant),] #dans le même ordre pour tous

    MSD_full <- MSD_full[,c("variant", "PX")]
    write.table(MSD_full, file = paste0("Traitwise/Pvalues_", TRAIT, ".csv"), sep=",", quote=F, row.names=F, col.names = T)

  #  print("start gzip")

    if(gzip == T){
      tryCatch(
        {        system(paste0("gzip Traitwise/Pvalues_", TRAIT, ".csv"))
        },
        error = function(cond) {
          message("gzip did not work but that's okay")
        })
    }



  }



  Analyse <- function(X){



    TRAIT = X
    print(paste0(TRAIT, " analysis"))

    path <- paste0("Traitwise/", list.files("Traitwise/", pattern = paste0("^Pvalues_",TRAIT,".csv")))

    MSD <- fread(path)

    TopVar <- na.omit(MSD$variant[MSD$PX < ThreshSelectionPvalues])
    TopVar_mem <- TopVar


    if(labelGWASsig == T){

    DepartMSD <- fread(paste0(sourceGWAS, list.files(sourceGWAS, pattern = paste0("^",TRAIT,"\\."))))
    if("pval" %in% colnames(DepartMSD)){
      Somme <- merge(MSD[, c("variant","PX")],  DepartMSD[, c("variant","pval")], by = "variant")
      Somme <- Somme[Somme$pval < ThreshSelectionPvalues | Somme$PX < ThreshSelectionPvalues,]
      TopVar <- Somme$variant
      TopVar_mem <- TopVar
      rm(Somme)
      rm(DepartMSD)
      gc()
    }else{  message("pval column absent from ", TRAIT ," GWAS input data. GWAS significant variants
                    that are not PRISM significant will not be labeled for pleiotropy.") }
    rm(DepartMSD)
    gc()
    }

    ### Vertical




    # We have to get all traits with a vertical effect on the trait of interest
    NTraits <- length(ListofTraits)
    TreshVpleio <- min(0.05,0.05/((NTraits * (NTraits-1) )/2))

    P1 <- ParametersTable[ParametersTable$X == TRAIT & ParametersTable$pval_ayx < TreshVpleio,]
    P2 <- ParametersTable[ParametersTable$Y == TRAIT & ParametersTable$pval_axy < TreshVpleio,]

    TraitsVPleio <-  c(P1$Y, P2$X)
    V_IDS <- NULL #If we don't enter the loop

    V_IDS_FULL <- c()




    for(VIND in TraitsVPleio){ #print(VIND)
      gc()
  #   print(VIND)



      if(file.exists(paste0("Pairwise/Likelihood_", TRAIT,"_", VIND,".csv.gz")) | file.exists(paste0("Pairwise/Likelihood_", TRAIT,"_", VIND,".csv"))){
        path <- paste0("Pairwise/", list.files("Pairwise/", pattern = paste0("Likelihood_", TRAIT,"_", VIND,".csv")))
        Omega_opti <- fread(path)
        Omega_opti <-  Omega_opti[Omega_opti$Om3 > Omega_opti$Om1 & Omega_opti$Om3 > Omega_opti$Om2 & Omega_opti$Om3 > Omega_opti$Om4 & Omega_opti$Om3 > Omega_opti$Om5 & Omega_opti$Om3 > Omega_opti$Om6 & Omega_opti$Om3 > Omega_opti$Om7,]
      }else{
        path <- paste0("Pairwise/", list.files("Pairwise/", pattern = paste0("Likelihood_", VIND,"_", TRAIT,".csv")))
        Omega_opti <- fread(path)
        Omega_opti <-  Omega_opti[Omega_opti$Om1 > Omega_opti$Om2 & Omega_opti$Om1 > Omega_opti$Om3 & Omega_opti$Om1 > Omega_opti$Om4 & Omega_opti$Om1 > Omega_opti$Om5 & Omega_opti$Om1 > Omega_opti$Om6 & Omega_opti$Om1 > Omega_opti$Om7,]
      }





      V_IDS <- data.frame( Omega_opti$INDEX)

      if(nrow(V_IDS)>0){V_IDS$Trait <- VIND
      row.names(V_IDS ) <- NULL

      V_IDS <- V_IDS[V_IDS$Omega_opti.INDEX %in% TopVar,]


      V_IDS_FULL <- rbind(V_IDS_FULL, V_IDS)}
    }



    quelVert <- function(VAR){
      Bop <- toString(V_IDS_FULL$Trait[which(V_IDS_FULL$Omega_opti.INDEX == VAR)])
      Bop <- str_replace_all(Bop , ", ", ":")
      Bop <- paste0("V:",Bop)
      return(Bop)
    }


    TopVar <- TopVar[TopVar %in% V_IDS_FULL$Omega_opti.INDEX]
    VertP <- mclapply(X = unique(TopVar), FUN = quelVert,  mc.cores = 1)
    VertP <- as.data.frame(unlist(VertP))
    VertP$variant <- unique(TopVar)
    if(nrow(VertP) > 0){colnames(VertP) <- c("VertP", "variant")}



    #### Effect on confounder prehension
    path <- paste0("Traitwise/", list.files("Traitwise/", pattern = paste0("Synthese_CatCharc_", TRAIT,".csv")))

    Conf <- fread(path, header = TRUE)



    MSD$Ori <- "No_supp_info"
    MSD$Ori2 <- MSD$Ori

    Conf <- Conf %>% filter_all(any_vars(. %in% c("Om2")))


    if(nrow(Conf) > 0){
      cols <- colnames(Conf)
      result <- data.frame(variant = Conf$variant,  Upleio = apply(Conf == "Om2", 1, function(x) toString(cols[x])))
      result$Upleio <- str_remove_all(result$Upleio, paste0(TRAIT,"_"))
      result$Upleio <- str_remove_all(result$Upleio, paste0("_",TRAIT))
      result$Upleio <- str_replace_all(result$Upleio, ", ", ":")
      result$Upleio <- paste0("U:",result$Upleio  )
      U_IDS <- Conf$variant




      MSD$Ori[match(U_IDS,MSD$variant)]  <-  result$Upleio
      MSD$Ori2[MSD$variant %in% U_IDS] <- "Confounder Pleiotropy"

    }



    MSD$Ori[match(VertP$variant,MSD$variant)]  <-  VertP$VertP
    #Vertical pleiotropy will be flagged as confounding pleiotropy (For another pair of traits, vertical is confounding)
    #So we input vertical last to overwrite confounding



    MSD$Ori2[MSD$variant %in% VertP$variant]  <- "Vertical Pleiotropy"
    MSD$Ori2[MSD$Ori == "No_supp_info"] <- "Direct Effect"
    MSD$Ori[MSD$Ori2 == "Direct Effect"] <- paste0("E:", TRAIT)
    MSD <- MSD[! is.na(MSD$PX),]
    MSD$Ori2[!MSD$variant %in%  TopVar_mem] <- "No Effect"
    MSD$Ori[!MSD$variant %in%  TopVar_mem] <- "No Effect"

    print(paste0(length(MSD$Ori2[MSD$Ori == "Direct Effect" & MSD$PX >= ThreshSelectionPvalues ] ), " variants are GWAS significant,
              but not PRISM significant, despite not being labeled for confounder or vertical pleiotropy. A high number could mean that PRISM
                 dos not accurately model trait ", TRAIT))
    MSD$Ori2[MSD$Ori == "Direct Effect" & MSD$PX >= ThreshSelectionPvalues ] <- "No Effect" #A variant cannot be direct if not PRISM sig
    MSD$Ori[MSD$Ori == "Direct Effect" & MSD$PX >= ThreshSelectionPvalues ] <- "No Effect" #Same


    MSD <- MSD[,c("variant", "PX", "Ori", "Ori2")]
    colnames(MSD) <- c("variant", "PvalPRISM", "FullPleio","SynthPleio")


    write.table(MSD, file = paste0("Results/Pleio_", TRAIT,".csv"), sep=",", quote=F, row.names=F, col.names = T)

  }



  mclapply(X = ListofTraits, FUN = SyntheseTraits, mc.cores = NbCores)




  mclapply(X = ListofTraits, FUN = Analyse, mc.cores = NbCores)

}
