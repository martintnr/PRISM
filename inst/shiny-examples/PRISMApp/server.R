library(ggplot2)
library(stats)
library(data.table)
library(parallel)
library(stringr)
library(rlist)
library(circlize)
library(splitstackshape)
library(visNetwork)
library(shiny)
library(shinydashboard)
library(RColorBrewer)
library(ggvis)
library(dplyr)

source("www/ggcircos_helpers.R")

# setwd('~/Partage/PleioMap/BioRes/Scripts/PRISM')

###########################################################################################################################################################################################
################################################################################# 0.1- Load PleioVar results ##############################################################################
###########################################################################################################################################################################################

if(! file.exists("www/VariantListPleio.Rdata")) {
    TopVariantsFiles <- list.files("www/Results_TopVar")
    VariantFileList <- lapply(TopVariantsFiles, function(file) system(paste0("tail -n +2 www/Results_TopVar/", file, " | awk -F \',\' \'{print $1}\'"), intern = TRUE))
    PleioType <- lapply(TopVariantsFiles, function(file) system(paste0("tail -n +2 www/Results_TopVar/", file, " | awk -F \',\' \'{print $4}\'"), intern = TRUE))

    VariantList <- cbind(SNPid = unlist(VariantFileList), PleioType = unlist(PleioType))
    VariantList <- by(unlist(VariantFileList), unlist(PleioType), table)
    VariantListPleio <- lapply(VariantList, function(variants) {
        VariantPleio <- cbind.data.frame(names(variants), t(simplify2array(strsplit(names(variants), split = ":"))), as.vector(variants))
        colnames(VariantPleio) = c("SNPid", "chromosome", "chromosome_end", "mutated_from_allele", "mutated_to_allele", "Number of pleiotropic assocations")
        VariantPleio$chromosome_start <- VariantPleio$chromosome_end <- as.integer(VariantPleio$chromosome_end)
        VariantPleio$chromosome <- as.character(VariantPleio$chromosome)
        return(VariantPleio)
    })
    names(VariantListPleio) <- names(VariantList)
    save(VariantListPleio, file = "www/VariantListPleio.Rdata")
} else {
    # List of trait files (61 elements)
    TopVariantsFiles <- list.files("www/Results_TopVar")
    # List of trait files, list of 61 elements with variant IDs
    VariantFileList <- lapply(TopVariantsFiles, function(file) system(paste0("tail -n +2 www/Results_TopVar/", file, " | awk -F \',\' \'{print $1}\'"), intern = TRUE))
    # Vect of all sig variants IDs
    AllVariants <- unique(unlist(VariantFileList))
    # List of the 3 pleiotropy types, with a table containing variant details for each pleio
    load("www/VariantListPleio.Rdata")
}

if(! file.exists("www/VariantListPleio_reduced.Rdata")){
    VariantFileList_reduced <- lapply(TopVariantsFiles, function(file) {
        cat(file)
        n <- as.numeric(gsub("(^ +)([0-9]+)( .*)", "\\2", system(paste0("wc -l www/Results_TopVar/", file), intern = TRUE)))
        system(paste0("tail -n +2 www/Results_TopVar/", file, " |  sort -k2 -g -t \',\' | head -n ", min(n-1, 500), " | awk -F \',\' \'{print $1}\'"), intern = TRUE)
    })
    VariantFileList_reduced <- unique(unlist(VariantFileList_reduced))
    length(unique(unlist(VariantFileList_reduced)))
    VariantListPleio_reduced <- lapply(VariantListPleio, function(dat) subset(dat, SNPid %in% VariantFileList_reduced))
    length(unique(unlist(lapply(VariantListPleio_reduced, `[`, 1))))
    save(VariantListPleio_reduced, file = "www/VariantListPleio_reduced.Rdata")
    
    # set.seed(123)
    # Percent <- .01
    # VariantListPleio_reduced <- lapply(VariantListPleio, function(dat) dat[sample(1:nrow(dat), size = Percent*nrow(dat)), ])
} else {
    load("www/VariantListPleio_reduced.Rdata")
}

###########################################################################################################################################################################################
################################################################################ 0.2- Build PleioVar Network ##############################################################################
###########################################################################################################################################################################################

getVariantNetwork <- function(Vex){

    getNetwork <- function(VariantPleioResult, PleioType){
        PleioResultSub <- subset(VariantPleioResult, SynthPleio == PleioType)
        
        # Direct effect
        if(PleioType == "No supplementary info" & any(PleioResultSub$SynthPleio == PleioType))
            Edges <- setnames(cbind.data.frame(PleioResultSub[, c("variant", "Trait1", "PvalPleioVar")], col = "#45709d"), c("from", "to", "value", "color"))
        
        # Vertical effect
        if(PleioType == "Suspected Vertical Pleiotropy" & any(VariantPleioResult$SynthPleio == PleioType)) {
            Edges <- rbind.data.frame(
                setnames(cbind.data.frame(PleioResultSub[, c("variant", "FullPleio", "PvalPleioVar")], col = "#82992a"), c("from", "to", "value", "color")),
                setnames(cbind.data.frame(PleioResultSub[, c("FullPleio", "Trait1", "PvalPleioVar")], col = "#82992a"), c("from", "to", "value", "color"))
            )
            # Remove vertical edges that are confirmed by direct effects
            EdgeVerif <- interaction(subset(VariantPleioResult, SynthPleio == "No supplementary info")[, c("variant", "Trait1")])
            Edges <- Edges[ ! interaction(Edges[, c("from", "to")]) %in% EdgeVerif, ]
        }
        # Network effect
        if(PleioType == "Detected Network Pleiotropy" & any(VariantPleioResult$SynthPleio == PleioType))
            Edges <- rbind.data.frame(
                setnames(cbind.data.frame(PleioResultSub[, c("variant", "PvalPleioVar")], to = paste("U", interaction(PleioResultSub[, c("Trait1", "FullPleio")], sep = "-"), sep = "-"), col = "#a23c33")[c(1, 3, 2, 4)], c("from", "to", "value", "color")),
                setnames(cbind.data.frame(PleioResultSub[, c("Trait1", "PvalPleioVar")], to = paste("U", interaction(PleioResultSub[, c("Trait1", "FullPleio")], sep = "-"), sep = "-"), col = "#a23c33")[c(3, 1, 2, 4)], c("from", "to", "value", "color")),
                setnames(cbind.data.frame(PleioResultSub[, c("FullPleio", "PvalPleioVar")], to = paste("U", interaction(PleioResultSub[, c("Trait1", "FullPleio")], sep = "-"), sep = "-"), col = "#a23c33")[c(3, 1, 2, 4)], c("from", "to", "value", "color"))
            )
        Edges$value <- -log10(as.numeric(Edges$value))
        return(Edges)
    }
    
    # Get variant results
    FullPleio <- do.call("rbind", lapply(TopVariantsFiles[grep(Vex, VariantFileList)], function(file){
        cmd <- system(paste0("grep ", Vex, " www/Results_TopVar/", file), intern = TRUE)
        if(length(cmd) > 0){
            cmd <- t(strsplit(cmd, split = ",")[[1]])
            colnames(cmd) <- c("variant", "PvalPleioVar", "FullPleio", "SynthPleio", "PvalGWAS")
            return(cbind.data.frame(Trait1 = gsub("^(Results_)(.+)(_.*\\.csv)$", "\\2", file), cmd))
        }
    }))
    FullPleio$FullPleio <- gsub("^[UV]:", "", FullPleio$FullPleio)
    FullPleio <- data.frame(cSplit(FullPleio, splitCols = "FullPleio", sep = ":", direction = "long", type.convert = TRUE))
    
    # Get Network
    PleioNetEdges <- do.call("rbind", lapply(unique(FullPleio$SynthPleio), getNetwork, VariantPleioResult = FullPleio))
    PleioNetEdges <- unique(PleioNetEdges)
    # PleioNetEdges <- PleioNetEdges[order(PleioNetEdges$value, decreasing = TRUE), ]
    # PleioNetEdges <- PleioNetEdges[! duplicated(paste(PleioNetEdges$from, PleioNetEdges$to)), ]
    
    # TraitDef <- fread("www/UKBB_TraitDef.txt", data.table = FALSE, select = c(1:3))
    Conf <- grep("U-", unique(unlist(PleioNetEdges[, c("from", "to")])), value = TRUE)
    TraitDefVex <- rbind.data.frame(TraitDef[, c(1, 3, 5, 2)], c(Category = "Variant", phenotype = Vex, description = Vex, color = "black"))
    PleioNetNodes <- setnames(subset(TraitDefVex, phenotype %in% unique(unlist(PleioNetEdges[, c("from", "to")]))), c("group", "id", "title", "color"))
    if(length(Conf) > 0){
        TraitDefVex <- rbind.data.frame(TraitDefVex, cbind.data.frame(Category = "Confounder", phenotype = Conf, description_short = "Confounder", Color = "#a23c33"))
        PleioNetNodes <- setnames(subset(TraitDefVex, phenotype %in% unique(unlist(PleioNetEdges[, c("from", "to")]))), c("group", "id", "title", "color"))
        PleioNetNodes$shape <- factor(as.numeric(PleioNetNodes$group == "Variant") + 10*as.numeric(PleioNetNodes$group == "Confounder"), labels = c("circle", "triangle", "box"))
        #PleioNetNodes$color <- factor(as.numeric(PleioNetNodes$group == "Variant") + 10*as.numeric(PleioNetNodes$group == "Confounder"), labels = c(NA, "black", "#a23c33"))
    } else{
        PleioNetNodes$shape <- factor(as.numeric(PleioNetNodes$group == "Variant"), labels = c("circle", "triangle"))
    }
    PleioNetNodes$label <- NA
    ########################################################################################################################################
    # deconvolution
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
            # paste0("Trait", i, " <- getEndTrait(PleioNetNodesIt = PleioNetNodes[! PleioNetNodes$group %in% c('Variant', 'Confounder') & ! PleioNetNodes$id %in% c(", paste0("Trait", 1:(i-1), collapse = ", "), "), ], PleioNetEdges) ; DeconvIndices <- c(DeconvIndices, which(PleioNetEdges$from %in% c(Vex, ", paste0("Trait", 1:(i-1), collapse = ", "), ") & PleioNetEdges$to %in% setdiff(PleioNetNodes$id[! PleioNetNodes$group %in% c('Variant', 'Confounder')], c(", paste0("Trait", 1:i, collapse = ", "), "))))")
            paste0("Trait", i, " <- getEndTrait(PleioNetNodesIt = PleioNetNodes[! PleioNetNodes$group %in% c('Variant', 'Confounder') & ! PleioNetNodes$id %in% c(", paste0("Trait", 1:(i-1), collapse = ", "), "), ], PleioNetEdges) ;", paste0(lapply(1:(i-1), function(j) paste0("Trait", j, "Sub <- PleioNetEdges[which(PleioNetEdges$from %in% Trait", j, " & PleioNetEdges$to %in% Trait", i, "), 'from']")), collapse = ";"), "; DeconvIndices <- c(DeconvIndices, which(PleioNetEdges$from %in% c(Vex, ", paste0("Trait", 1:(i-1), "Sub", collapse = ", "), ") & PleioNetEdges$to %in% setdiff(PleioNetNodes$id[! PleioNetNodes$group %in% c('Variant', 'Confounder')], c(", paste0("Trait", 1:i, collapse = ", "), "))))")
        )
        
        ########## Principle
            # Get end trait at step 3 for instance
        # Trait3 <- getEndTrait(PleioNetNodesIt = PleioNetNodes[! PleioNetNodes$group %in% c('Variant', 'Confounder') & ! PleioNetNodes$id %in% c(Trait1, Trait2), ], PleioNetEdges)
            # Reassess traits sets at previous steps
        # Trait2Sub <- PleioNetEdges[which(PleioNetEdges$from %in% Trait2 & PleioNetEdges$to %in% Trait3), "from"]
        # Trait1Sub <- PleioNetEdges[which(PleioNetEdges$from %in% Trait1 & PleioNetEdges$to %in% Trait3), "from"] # Check that all traits from trait 1 go to trait 2
            # Remove edges
        # DeconvIndices3 <- which(PleioNetEdges$from %in% c(Vex, Trait1Sub, Trait2Sub) & PleioNetEdges$to %in% setdiff(PleioNetNodes$id[! PleioNetNodes$group %in% c('Variant', 'Confounder')], c(Trait1, Trait2, Trait3)))
        ######################
        
        eval(parse(text = paste(DeconvProcess, collapse = ";")))
        DeconvIndices <- unique(DeconvIndices)
        PleioNetEdgesDeconv <- PleioNetEdges[-DeconvIndices, ]
        
        # Remove duplicated edges
        PleioNetEdgesDeconv <- PleioNetEdgesDeconv[order(PleioNetEdgesDeconv$value, decreasing = TRUE), ]
        PleioNetEdgesDeconv <- PleioNetEdgesDeconv[! duplicated(paste(PleioNetEdgesDeconv$from, PleioNetEdgesDeconv$to)), ]
    } else {
        PleioNetEdgesDeconv <- PleioNetEdges
    }
    return(list(Edges = PleioNetEdges, EdgesDeconv = PleioNetEdgesDeconv, Nodes = PleioNetNodes))
}

###########################################################################################################################################################################################
################################################################################ 0.3- Load trait information ##############################################################################
###########################################################################################################################################################################################

# Trait relationship
TraitCharacteristics <- fread("www/Table_Parameters_alltraits_pval.csv", data.table = FALSE)
TraitDef <- fread("www/UKBB_TraitDef.txt", data.table = FALSE, select = 1:5)
TraitCharacteristics <- merge(TraitDef[, c(1,2,3,5)], TraitCharacteristics, by.x = "phenotype", by.y = "Y")
colnames(TraitCharacteristics)[1:4] <- paste(c(colnames(TraitCharacteristics)[1:3], "description"), "Y", sep = "_")
TraitCharacteristics <- merge(TraitDef[, c(1,2,3,5)], TraitCharacteristics, by.x = "phenotype", by.y = "X")
colnames(TraitCharacteristics)[1:4] <- paste(c(colnames(TraitCharacteristics)[1:3], "description"), "X", sep = "_")
TraitCharacteristics$description_X <- factor(TraitCharacteristics$description_X, levels = unique(TraitCharacteristics$description_X[order(TraitCharacteristics$Category_X)]))
TraitCharacteristics$description_Y <- factor(TraitCharacteristics$description_Y, levels = unique(TraitCharacteristics$description_Y[order(TraitCharacteristics$Category_Y)]))
nbTraits <- length(unique(TraitCharacteristics$phenotype_Y))
grid.col <- TraitCharacteristics$Color_X
names(grid.col) <- TraitCharacteristics$description_X

###########################################################################################################################################################################################
######################################################################### 0.4- Load chromosome representation data ########################################################################
###########################################################################################################################################################################################

PleioColors <- c("#45709d", "#82992a", "#a23c33")
names(PleioColors) <- c("No supplementary info", "Suspected Vertical Pleiotropy", "Detected Network Pleiotropy")

scalingFactor <- 10

chroms <- as.character(c(1:22))
lengths <- c(
    249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
    159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
    115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
    59128983, 63025520, 48129895, 51304566
)

radians_f <- create_radians(chroms[1:21], lengths[1:21])
radians_m <- create_radians(chroms, lengths)

# chrom <- 21
# VariantListNR <- sort(grep(paste0("^", chrom, "\\:"), as.character(unlist(lapply(VariantListPleio, `[[`, 1))), value = TRUE))

###########################################################################################################################################################################################
################################################################################### 1- Start shiny app ####################################################################################
###########################################################################################################################################################################################

# showReactLog()

shinyServer(function(input, output, session){
    
    
    #######################################################################################################################################################################################
    ###################################################################################### About tab ######################################################################################
    #######################################################################################################################################################################################
    
    output$TournaireEtAl <- renderUI({
        tagList(a(
            "Tournaire et al. 2024.",
            # "Inferring genetic variant causal network by leveraging pleiotropy",
            href = "https://www.medrxiv.org/content/10.1101/2024.06.01.24308193v1",
            target = "_blank"
        ))
    })
    output$github <- renderUI({
        tagList(a(
            "martintnr/PRISM", 
            href = "https://github.com/martintnr/PleioVar",
            target = "_blank"
        ))
    })
    output$email <- renderUI({
        tagList(a(
            h5("marie.verbanck [at] u-paris.fr"), 
            href = "mailto:marie.verbanck@u-paris.fr")
        )
    })
    output$info <- renderUI({
        tagList(a(
            h5("Further information on the GEnetics Machine Learning lab"),
            href = "http://marie.verbanck.free.fr/",
            target = "_blank")
        )
    })
    
    output$AboutPage <- renderUI({
        tagList(
            a("Back to the about tab", onclick = "openTab('About')"),
            tags$script(HTML("
				var openTab = function(tabName){
		  			$('a', $('.sidebar')).each(function() {
						if(this.getAttribute('data-value') == tabName) {
							this.click()
						};
					});
				}
			")
            ))
    })

    #######################################################################################################################################################################################
    ############################################################################ A- Trait Causal relationships ############################################################################
    #######################################################################################################################################################################################
    
    output$TraitCausal <- renderPlot({
        
        circos.clear()
        circos.par(gap.after = rep(0.5, nbTraits))
        chordDiagram(
            cbind.data.frame(TraitCharacteristics[, c("description_X", "description_Y")], replace(TraitCharacteristics$axy, TraitCharacteristics$pval_axy > as.numeric(input$PvalCutoff), 0)), 
            order = unique(TraitCharacteristics$description_X[order(TraitCharacteristics$Category_X)]),
            grid.col = grid.col,
            directional = 1,
            direction.type = "arrows",link.arr.type = "big.arrow", link.arr.length = 0.2,
            annotationTrack = c("grid"), annotationTrackHeight = c(0.35,0.1),
            reduce = 0
        )
        
        circos.track(track.index = 1, panel.fun = function(x, y) {
            xlim = get.cell.meta.data("xlim")
            ylim = get.cell.meta.data("ylim")
            sector.index = get.cell.meta.data("sector.index")
            circos.text(mean(xlim), mean(ylim), sector.index, col = "white", cex = 0.9, facing = "clockwise", niceFacing = TRUE)
        }, bg.border = NA)
    })
    
    #######################################################################################################################################################################################
    ############################################################################# B- Chromosome Representation ############################################################################
    #######################################################################################################################################################################################
    re_values <- reactiveValues(
        scaling_factors = rep(1, 22),
        previous_radians = radians_f,
        chrom_clicked = "1", 
        chroms_selected = NULL
    )

    VariantList <- reactive({
        if(is.null(re_values$chroms_selected))
            return(NULL)
        VariantList <- sort(grep(paste0("^", re_values$chroms_selected, "\\:"), as.character(unlist(lapply(VariantListPleio_reduced, `[[`, 1))), value = TRUE))
    })
    
    radians <- reactive({
        rads <- create_radians(chroms, lengths * re_values$scaling_factors)
        isolate(mid <- mean(rads[names(rads) == re_values$chrom_clicked]))
        isolate(prev_mid <- mean(re_values$previous_radians[names(rads) == re_values$chrom_clicked]))
        offset <- mid - prev_mid
        rads - offset
    })
    
    track_radians <- reactive({
        create_track_radians(radians(), points_per_track = rep(20, length(radians())))
    })
    
    seq_df <- reactive({
        scale <- re_values$scaling_factors
        create_seq_df(radians(), scale = scale)
    })
    
    snp_plot_data <- reactive({
        snp <- VariantListPleio_reduced[input$PleioType][[1]]
        snp <- snp[snp$`Number of pleiotropic assocations` > 1, ]
        
        snp_data <- snp %>%
            group_by(
                SNPid, chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele
            ) %>%
            summarise(transcripts = `Number of pleiotropic assocations`)
        snp_data$chromosome <- as.character(snp_data$chromosome)
        points <- fit_points(snp_data$chromosome, snp_data$chromosome_start, snp_data$transcripts, 0.8, 0.6, seq_df(), metadata = snp_data[, c("SNPid", "transcripts", "chromosome","chromosome_start", "chromosome_end", "mutated_from_allele", "mutated_to_allele")], min_value = 0, max_value = max(snp_data$transcripts) + 1
        )
        points$id <- paste0("snp", 1:nrow(points))
        points
    })
    
    text_df <- reactive({
        seq_df() %>% mutate(theta = (seq_start + seq_end) / 2, r = 1.05)
    })
    
    tooltip_fun <- function(data, session) {
        if ("id" %in% names(data)) {
            tt_data <- snp_plot_data()
            row <- tt_data[tt_data$id == data$id, ]
            paste0(
                "<b>Variant: ", row$SNPid, "</b><br>",
                "Position: ", row$chromosome_start, "<br>",
                "Base Change: ", row$mutated_from_allele, ">", row$mutated_to_allele, "<br>",
                "Number of pleiotropic associations: ", row$transcripts, "<br>"
            )
        }
    }
    
    # tooltip_click_fun <- function(data) {
    #   str(data)
    # }
    
    click_handle <- function(data, location, session) {
        if (is.null(data)) {
            return(NULL)
        }
        
        isolate(re_values$chrom_clicked <- data$group)
        isolate(re_values$previous_radians <- radians())
        isolate(re_values$scaling_factors[which(chroms == data$group)] <- ifelse(re_values$scaling_factors[which(chroms == data$group)] == 1, scalingFactor, 1))
        isolate(re_values$chroms_selected <- data$group)
        print(data)
    }
    
    fill_range <- stroke_range <- c(
        # chromosome colours from Circos
        "#996600", "#666600", "#99991E", "#CC0000", "#FF0000", "#FF00CC", "#FFCCCC", "#FF9900", "#FFCC00",
        "#FFFF00", "#CCFF00", "#00FF00", "#358000", "#0000CC", "#6699FF", "#99CCFF", "#00FFFF", "#CCFFFF",
        "#9900CC", "#CC33FF", "#CC99FF", "#666666"
    )
    
    fill_domain <- stroke_domain <- c(1:22)
    
    add_tooltip <- function(vis, html, on = c("hover", "click")) {
        on <- match.arg(on)
        
        show_tooltip2 <- function(data, location, session, ...) {
            if (is.null(data)) {
                hide_tooltip(session)
                return()
            }
            
            html <- html(data)
            if (is.null(html)) {
                hide_tooltip(session)
            } else {
                show_tooltip(session, location$x + 5, location$y + 5, html)
            }
        }
        hide_tooltip2 <- function(session) {
            hide_tooltip(session)
        }
        
        switch(on,
               click = handle_click(vis, show_tooltip2),
               hover = handle_hover(vis, show_tooltip2)
        )
    }
    
    ggvis() %>%
        add_track(track_radians, 1, 0.9, fill = ~group, stroke = ~group, fillOpacity := 0.7, fillOpacity.hover := 1) %>%
        add_track(track_radians, 0.8, 0.6, strokeOpacity := 0.5, stroke := "black", strokeWidth := 0.5) %>%
        add_circles(track_radians, seq(0.6, 0.8, length.out = 7)[-c(1, 7)], strokeOpacity := 0.3, strokeWidth := 0.5) %>%
        # Layer with points
        layer_points(
            data = snp_plot_data, ~ sin(theta) * r, ~ cos(theta) * r, size := 10,
            key := ~id, 
            size.hover := 30, strokeOpacity.hover := 1, fill := reactive({as.character(PleioColors[input$PleioType])}),
            stroke := reactive({as.character(PleioColors[input$PleioType])}), strokeWidth := 0.5
        ) %>%
        layer_text(
            data = text_df, ~ sin(theta) * r, ~ cos(theta) * r, text := ~seq, align := "center", baseline := "middle",
            angle := ~ 180 * (theta - pi * (cos(theta) < 0)) / pi
        ) %>%
        add_tooltip(tooltip_fun, "hover") %>%
        # add_tooltip(tooltip_click_fun, "click") %>%
        handle_click(click_handle) %>% # Allow to unfold the chromosome
        scale_ordinal("fill", domain = fill_domain, range = fill_range) %>%
        scale_ordinal("stroke", domain = stroke_domain, range = stroke_range) %>%
        scale_nominal("shape", domain = c("Yes", "No"), range = c("cross", "circle")) %>%
        scale_nominal("size", domain = c("Yes", "No"), range = c(20, 10)) %>%
        scale_nominal("opacity", domain = c("Yes", "No"), range = c(1, 0)) %>%
        hide_axis("x") %>%
        hide_axis("y") %>%
        hide_legend(c("fill", "stroke", "shape", "size")) %>%
        set_options(hover_duration = 0, width = 1100, height = 1100, keep_aspect = TRUE, duration = 1000) %>%
        bind_shiny("ChromosomePlot")

    #######################################################################################################################################################################################
    ################################################################################### C- Variant network ################################################################################
    #######################################################################################################################################################################################
    
    observeEvent(input$GoToFullResults, {
        updateTabItems(session, "tabs", "VariantNetFull")
    })
    
    output$VariantSelect <- renderUI({
        VariantList <- VariantList()
        if(is.null(VariantList))
            return(NULL)
        box(
            title = "Explore the network visualization for the selected genetic variant",
            selectizeInput("variant", label = "", choices = VariantList), # Set choices to null so that it load more quickly and then update selectizeInput in server with the choices
            status = "primary",
            solidHeader = TRUE, 
            width = 8
        )
    })
    output$NbVariants <- renderUI({
        VariantList <- VariantList()
        if(is.null(VariantList))
            return(NULL)
        valueBox(length(VariantList), paste("genetic variants in chromosome", re_values$chroms_selected), icon = icon("dna"), width = 4)
    })
    
    output$VariantNetwork <- renderVisNetwork({
        Net <- getVariantNetwork(input$variant)
        if(is.null(Net))
            return(NULL)
        visNetwork(Net$Nodes, Net$Edges, width = 300, height = 600) %>% 
            visEdges(arrows = 'to', scaling = list(min = 2, max = 2)) %>%
            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group")  %>%
            visLegend(useGroups = FALSE, addNodes = cbind.data.frame(setnames(unique(Net$Nodes[, c("group", "shape", "color")]), c("label", "shape", "color")), font.color = "lightgrey"))
    })
    
    output$VariantNetworkDeconv <- renderVisNetwork({
        Net <- getVariantNetwork(input$variant)
        if(is.null(Net))
            return(NULL)
        visNetwork(Net$Nodes, Net$EdgesDeconv, width = 300, height = 600) %>% 
            visEdges(arrows = 'to', scaling = list(min = 2, max = 2)) %>%
            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
            visLegend(useGroups = FALSE, addNodes = cbind.data.frame(setnames(unique(Net$Nodes[, c("group", "shape", "color")]), c("label", "shape", "color")), font.color = "lightgrey"))
    })
    
    
    output$VariantNetworkTabs <- renderUI({
        if(is.null(re_values$chroms_selected))
            return(NULL)
        box(
            tabsetPanel(type = "tabs",
                        tabPanel("Deconvoluted",
                                 visNetworkOutput("VariantNetworkDeconv")
                        ),
                        tabPanel("Raw network",
                                 visNetworkOutput("VariantNetwork"),
                        )
            ),
            title = "Network visualization for the selected variant",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12
        )
    })

    #######################################################  
    output$VariantNetworkFull <- renderVisNetwork({
        if(! input$ManualVariant %in% AllVariants)
            return(NULL)
        Net <- getVariantNetwork(input$ManualVariant)
        visNetwork(Net$Nodes, Net$Edges, width = 300, height = 600) %>% 
            visEdges(arrows = 'to', scaling = list(min = 2, max = 2)) %>%
            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group")  %>%
            visLegend(useGroups = FALSE, addNodes = cbind.data.frame(setnames(unique(Net$Nodes[, c("group", "shape", "color")]), c("label", "shape", "color")), font.color = "lightgrey"))
    })
    
    output$VariantNetworkDeconvFull <- renderVisNetwork({
        if(! input$ManualVariant %in% AllVariants)
            return(NULL)
        Net <- getVariantNetwork(input$ManualVariant)
        visNetwork(Net$Nodes, Net$EdgesDeconv, width = 300, height = 600) %>% 
            visEdges(arrows = 'to', scaling = list(min = 2, max = 2)) %>%
            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
            visLegend(useGroups = FALSE, addNodes = cbind.data.frame(setnames(unique(Net$Nodes[, c("group", "shape", "color")]), c("label", "shape", "color")), font.color = "lightgrey"))
    })
    
    output$VariantNetworkTabsFull <- renderUI({
        if(! input$ManualVariant %in% AllVariants){
            infoBox(
                title = "Warning",
                "The variant has not been found in our database.",
                color = "red",
                width = 12,
                icon = icon("triangle-exclamation"),
                fill = TRUE
            )
        } else {
            box(
                tabsetPanel(type = "tabs",
                            tabPanel("Deconvoluted",
                                     visNetworkOutput("VariantNetworkDeconvFull")
                            ),
                            tabPanel("Raw network",
                                     visNetworkOutput("VariantNetworkFull"),
                            )
                ),
                title = "Network visualization for the selected variant",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 12
            )
        }
    })
    
})
