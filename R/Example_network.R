#' Graph network from the example data
#'
#' `Example_network()` returns a ggplot of the results obtained from the example.
#'
#' @param Variant The variant of interest.
#' @param ListofTraits All traits in the analysis (necessary to calculate the p-value threshold)
#' @param ThreshSelectionPvalues P-value threshold to select PRISM top variants, if you do not want
#' the default Bonferroni correction.
#'@param ParametersTable The parameters table for each pair of traits.
#'
#' @return a deconvoluted causal variant network
#' @export
#'
#' @examples
Example_network <- function(ListofTraits, Variant, ParametersTable, ThreshSelectionPvalues = 5e-08/(length(ListofTraits)-1)){




  TraitDef <- data.frame(Category = c("A","B","B", "B", "B", "C","C", "C", "C", "D",
                                      "E","E","E","E","E","E","E","E") ,
                         phenotype = c("A", #Anciens traits triés
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
                                         "E6",
                                         "E7",
                                         "E8"))

  TraitDef$description <- TraitDef$phenotype


  message("Aggregating top variants...")


  if(!file.exists("Networks/")){system("mkdir Networks")}


  VariantFileList <- lapply(ListofTraits, function(TRAIT){
    if(file.exists( paste0("Networks/TopVariants_", TRAIT,".csv"))){
      print(paste0("Networks/TopVariants_", TRAIT,".csv", " already exists."))
      TAB <- fread(paste0("Networks/TopVariants_", TRAIT,".csv"))
      return(as.character(TAB$variant))


    }else{


      TAB <- fread(paste0("Results/Pleio_", TRAIT,".csv"))
      TAB <- TAB[TAB$PvalPRISM < ThreshSelectionPvalues,]
      write.table(TAB, file = paste0("Networks/TopVariants_", TRAIT,".csv"), sep = ",", row.names = F)
      return(as.character(TAB$variant))

    }
  })

  message("Done")

  TopVariantsFiles <- lapply(ListofTraits, function(TRAIT){
    return(paste0("TopVariants_", TRAIT,".csv"))
  })



  Vex <- Variant

  GetAnyNetworks <- function(Vex){

    getNetwork <- function(VariantPleioResult, PleioType){

      PleioResultSub <- subset(VariantPleioResult, SynthPleio == PleioType)

      # Direct effect
      if(PleioType == "No supplementary info" & any(PleioResultSub$SynthPleio == PleioType))
        Edges <- setnames(cbind.data.frame(PleioResultSub[, c("variant", "Trait1", "PvalPRISM")], col = "#45709d"), c("from", "to", "value", "color"))

      # Vertical effect
      if(PleioType == "Suspected Vertical Pleiotropy" & any(VariantPleioResult$SynthPleio == PleioType)) {
        Edges <- rbind.data.frame(
          setnames(cbind.data.frame(PleioResultSub[, c("variant", "FullPleio", "PvalPRISM")], col = "#82992a"), c("from", "to", "value", "color")),
          setnames(cbind.data.frame(PleioResultSub[, c("FullPleio", "Trait1", "PvalPRISM")], col = "#82992a"), c("from", "to", "value", "color"))
        )
        # Remove vertical edges that are confirmed by direct effects
        EdgeVerif <- interaction(subset(VariantPleioResult, SynthPleio == "No supplementary info")[, c("variant", "Trait1")])
        Edges <- Edges[ ! interaction(Edges[, c("from", "to")]) %in% EdgeVerif, ]
      }
      # Network effect
      if(PleioType == "Detected Network Pleiotropy" & any(VariantPleioResult$SynthPleio == PleioType))
        Edges <- rbind.data.frame(
          setnames(cbind.data.frame(PleioResultSub[, c("variant", "PvalPRISM")], to = paste("U", interaction(PleioResultSub[, c("Trait1", "FullPleio")], sep = "-"), sep = "-"), col = "#a23c33")[c(1, 3, 2, 4)], c("from", "to", "value", "color")),
          setnames(cbind.data.frame(PleioResultSub[, c("Trait1", "PvalPRISM")], to = paste("U", interaction(PleioResultSub[, c("Trait1", "FullPleio")], sep = "-"), sep = "-"), col = "#a23c33")[c(3, 1, 2, 4)], c("from", "to", "value", "color")),
          setnames(cbind.data.frame(PleioResultSub[, c("FullPleio", "PvalPRISM")], to = paste("U", interaction(PleioResultSub[, c("Trait1", "FullPleio")], sep = "-"), sep = "-"), col = "#a23c33")[c(3, 1, 2, 4)], c("from", "to", "value", "color"))
        )

      Edges$value <- -log10(as.numeric(Edges$value))
      return(Edges)
    }

    Quelfichier <- unlist(lapply(X = c(1:length(VariantFileList)), FUN = function(number){
      if(Vex %in% VariantFileList[[number]]){return(number)}
    }) )

    FullPleio <- do.call("rbind", lapply(TopVariantsFiles[Quelfichier], function(file){
      cmd <- system(paste0("grep ", Vex, " Networks/", file), intern = TRUE)
      if(length(cmd) > 0){
        cmd <- t(strsplit(cmd, split = ",")[[1]])
        colnames(cmd) <- c("variant", "PvalPRISM", "FullPleio", "SynthPleio")
        return(cbind.data.frame(Trait1 = gsub("^(TopVariants_)(.+)(_.*\\.csv)$", "\\2", file), cmd))
      }
    }))

    FullPleio$FullPleio <- gsub("\"", "", FullPleio$FullPleio)
    FullPleio$SynthPleio <- gsub("\"", "", FullPleio$SynthPleio)

    FullPleio$Trait1 <- gsub("TopVariants_", "", FullPleio$Trait1)
    FullPleio$Trait1 <- gsub(".csv", "",  FullPleio$Trait1)

    FullPleio$variant <- gsub("\"", "",  FullPleio$variant)



    FullPleio$FullPleio <- gsub("^[UV]:", "", FullPleio$FullPleio)
    FullPleio <- data.frame(cSplit(FullPleio, splitCols = "FullPleio", sep = ":", direction = "long", type.convert = TRUE))
    FullPleio <- FullPleio[order(as.numeric(FullPleio$PvalPRISM), decreasing = T), ]

    # Remove duplicated confounders
    toremove <- c()
    Conf1 <- c()
    Conf2 <- c()
    for(L in c(1:nrow(FullPleio))){

      if(FullPleio$SynthPleio[L] == "Confounder Pleiotropy"){
        if((FullPleio$Trait1[L] %in% Conf1 &  FullPleio$FullPleio[L]%in% Conf2) | (FullPleio$Trait1[L] %in% Conf2 &  FullPleio$FullPleio[L]%in% Conf1)){
          toremove <- c(toremove, L)}else{
            Conf1 <- c(Conf1, FullPleio$Trait1[L])
            Conf2 <- c(Conf2, FullPleio$FullPleio[L]) }
      }
    }
    if(!is.null(toremove)){FullPleio <- FullPleio[-toremove,]}



    # Get Network
    PleioNetEdges <- do.call("rbind", lapply(unique(FullPleio$SynthPleio), getNetwork, VariantPleioResult = FullPleio))
    PleioNetEdges <- unique(PleioNetEdges)
    Conf <- grep("U-", unique(unlist(PleioNetEdges[, c("from", "to")])), value = TRUE)
    TraitDef <- rbind.data.frame(TraitDef, c(Category = "Variant", phenotype = Vex, description = Vex))
    PleioNetNodes <- setnames(subset(TraitDef, phenotype %in% unique(unlist(PleioNetEdges[, c("from", "to")]))), c("group", "id", "title"))
    PleioNetNodes$label <- NA
    if(length(Conf) > 0){
      TraitDef <- rbind.data.frame(TraitDef, cbind.data.frame(Category = "Confounder", phenotype = Conf, description = "Confounder"))
      PleioNetNodes <- setnames(subset(TraitDef, phenotype %in% unique(unlist(PleioNetEdges[, c("from", "to")]))), c("group", "id", "title"))
      PleioNetNodes$shape <- factor(as.numeric(PleioNetNodes$group == "Variant") + 10*as.numeric(PleioNetNodes$group == "Confounder"), labels = c("circle", "triangle", "box"))
      PleioNetNodes$color <- factor(as.numeric(PleioNetNodes$group == "Variant") + 10*as.numeric(PleioNetNodes$group == "Confounder"), labels = c(NA, "black", "#a23c33"))
    } else{
      PleioNetNodes$shape <- factor(as.numeric(PleioNetNodes$group == "Variant"), labels = c("circle", "triangle"))
      PleioNetNodes$color <- factor(as.numeric(PleioNetNodes$group == "Variant"), labels = c(NA, "black"))
    }


    server <- function(input, output) {
      output$VariantNetwork <- renderVisNetwork({
        visNetwork(PleioNetNodes, PleioNetEdges, width = 300, height = 600) %>%
          visEdges(arrows = 'to', scaling = list(min = 2, max = 2)) %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group")
        # %>% visLegend()
      })
    }

    ui <- fluidPage(
      visNetworkOutput("VariantNetwork")
    )

    #shinyApp(ui = ui, server = server)

    ########################################################################################################################################

    # deconv
    # - if vertical effect = need to check for direct effect of vertical trait
    # - if multiple vertical effects: need to cut conditional effects

    getEndTrait <- function(PleioNetNodesIt, PleioNetEdges){
      unlist(lapply(PleioNetNodesIt$id, function(trait){
        if(nrow(PleioNetEdges[PleioNetEdges$from %in% PleioNetNodesIt$id[! PleioNetNodesIt$group %in% c("Variant", "Confounder")] & PleioNetEdges$to == trait, ]) == 0)
          return(trait)
      }))
    }

    # Step-by-step procedure to remove edges from the variant to the end traits
    Trait1 <- getEndTrait(PleioNetNodesIt = PleioNetNodes[! PleioNetNodes$group %in% c("Variant", "Confounder"), ], PleioNetEdges)
    DeconvIndices <- which(PleioNetEdges$from == Vex & PleioNetEdges$to %in% setdiff(PleioNetNodes$id[! PleioNetNodes$group %in% c("Variant", "Confounder")], Trait1))
    if(sum(DeconvIndices) > 0) {
      DeconvProcess <- lapply(2:(length(PleioNetNodes$id) -2), function(i)
        paste0("Trait", i, " <- getEndTrait(PleioNetNodesIt = PleioNetNodes[! PleioNetNodes$group %in% c('Variant', 'Confounder') & ! PleioNetNodes$id %in% c(", paste0("Trait", 1:(i-1), collapse = ", "), "), ], PleioNetEdges) ; DeconvIndices <- c(DeconvIndices, which(PleioNetEdges$from %in% c(Vex, ", paste0("Trait", 1:(i-1), collapse = ", "), ") & PleioNetEdges$to %in% setdiff(PleioNetNodes$id[! PleioNetNodes$group %in% c('Variant', 'Confounder')], c(", paste0("Trait", 1:i, collapse = ", "), "))))")
      )
      eval(parse(text = paste(DeconvProcess, collapse = ";")))
      DeconvIndices <- unique(DeconvIndices)
      PleioNetEdgesDeconv <- PleioNetEdges[-DeconvIndices, ]

      # Remove duplicated edges
      PleioNetEdgesDeconv <- PleioNetEdgesDeconv[order(PleioNetEdgesDeconv$value, decreasing = TRUE), ]
      PleioNetEdgesDeconv <- PleioNetEdgesDeconv[! duplicated(paste(PleioNetEdgesDeconv$from, PleioNetEdgesDeconv$to)), ]
    } else {
      PleioNetEdgesDeconv <- PleioNetEdges
    }



    # Shiny representation
    server <- function(input, output) {
      output$VariantNetwork <- renderVisNetwork({
        visNetwork(PleioNetNodes, PleioNetEdgesDeconv, width = 300, height = 600) %>%
          visEdges(arrows = 'to', scaling = list(min = 2, max = 2)) %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group")
        # %>% visLegend()
      })
    }

    ui <- fluidPage(
      visNetworkOutput("VariantNetwork")
    )


   print( shinyApp(ui = ui, server = server) )

    return(PleioNetEdgesDeconv)
  }


  Network_asked <- GetAnyNetworks(Vex)

  return(Network_asked)
}

