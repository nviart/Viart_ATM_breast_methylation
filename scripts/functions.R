#seed = seed

# Retrieve group of annotations
getGroups <- function(annotations){
    inactive <- c("T0072", "T0249", "T0076", "T0077R", "T0078",
                  "T0099", "T0220", "T0248", "T0009", "T0247",
                  "M5137110001", "M5769110001", "S2006101", "T0003", "T0141", "T0181")
    
    annotations[["Inactive"]] <- ifelse(annotations$Sample_Name %in% inactive & annotations$Sample_Group == "ATM",
                                        "inactive ATM",
                                 ifelse((! annotations$Sample_Name %in% inactive) & annotations$Sample_Group == "ATM", "active ATM",
                                        "Control")
                                 )
    
    #l1 <- c("M5137110001", "M5616110001", "M5738110001", "M5769110001", "S2006101", "S2013401", "S2022501", "S2029601", "S5033210001")
    #l2 <- c("T0544", "T0545", "T0546", "T0547", "T0548", "T0549", "T0550", "T0551", "T0552")
    
    #for (i in 1:length(l1)) {
    #    print(l1[i])
    #    test <- data.frame(lapply(test, function(x) {
    #        gsub(l1[i], l2[i], x)
    #    }))
    #}    
    
    samples_Curie <- subset(annotations, annotations$Sample_Group == "ATM" & annotations$ER_status == "1")
    samples_control <- subset(annotations, annotations$Sample_Group == "Non.ATM" & annotations$ER_status == "1")
    annotations_A0A2 <- rbind(samples_Curie, samples_control)


    samples_Curie <- subset(annotations, annotations$Sample_Group == "ATM" & annotations$ClinVar %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"))
    samples_control <- subset(annotations, annotations$Sample_Group == "Non.ATM")
    annotations_A1 <- rbind(samples_Curie, samples_control)
    
    samples_Curie <- subset(annotations, annotations$Sample_Group == "ATM" & annotations$ER_status == "1" & annotations$ClinVar %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"))
    samples_control <- subset(annotations, annotations$Sample_Group == "Non.ATM" & annotations$ER_status == "1")
    annotations_A2 <- rbind(samples_Curie, samples_control)
    
    samples_Curie <- subset(annotations, annotations$Sample_Group == "ATM" & annotations$subtype %in% c("Luminal A", "Luminal B", "Luminal A/B") & annotations$ClinVar %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"))
    samples_control <- subset(annotations, annotations$Sample_Group == "Non.ATM" & annotations$subtype %in% c("Luminal A", "Luminal B", "Luminal A/B"))
    annotations_A3 <- rbind(samples_Curie, samples_control)
    
    samples_Curie <- subset(annotations, annotations$Sample_Group == "ATM" & annotations$subtype == "Luminal B" & annotations$ClinVar %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"))
    samples_control <- subset(annotations, annotations$Sample_Group == "Non.ATM" & annotations$subtype == "Luminal B")
    annotations_A4 <- rbind(samples_Curie, samples_control)

    samples_Curie <- subset(annotations, annotations$Sample_Group == "ATM" & annotations$ClinVar %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"))
    samples_control <- subset(annotations, annotations$Sample_Group == "Non.ATM" & annotations$Diagnosis_age_categorical != "Young")
    annotations_AOld <- rbind(samples_Curie, samples_control)

    samples_Curie <- subset(annotations, annotations$Sample_Group == "ATM" & annotations$ClinVar %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic") & annotations$ER_status == "1")
    samples_control <- subset(annotations, annotations$Sample_Group == "Non.ATM" & annotations$Diagnosis_age_categorical != "Young" & annotations$ER_status == "1")
    annotations_AOldER <- rbind(samples_Curie, samples_control)

    samples_Curie <- subset(annotations, annotations$Sample_Group == "ATM" & ! annotations$ClinVar %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic") )
    samples_control <- subset(annotations, annotations$Sample_Group == "Non.ATM")
    annotations_AVUS <- rbind(samples_Curie, samples_control)

    samples_Curie <- subset(annotations, annotations$Sample_Group == "ATM" & ! annotations$ClinVar %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic") & annotations$ER_status == "1")
    samples_control <- subset(annotations, annotations$Sample_Group == "Non.ATM" & annotations$ER_status == "1")
    annotations_AVUSA2 <- rbind(samples_Curie, samples_control)
    
    annotations_A2Inactive <- annotations_A2[(annotations_A2$Sample_Group == "ATM" & annotations_A2$Inactive == "inactive ATM") |
                                             annotations_A2$Sample_Group == "Non.ATM" , ]

    return(list(
        "annotations" = annotations,
        "annotations_A0A2" = annotations_A0A2,
        "annotations_A1" = annotations_A1,
        "annotations_A2" = annotations_A2,
        "annotations_A3" = annotations_A3,
        "annotations_A4" = annotations_A4,
        "annotations_AOld" = annotations_AOld,
        "annotations_AOldER" = annotations_AOldER,
        "annotations_AVUS" = annotations_AVUS,
        "annotations_AVUSA2" = annotations_AVUSA2,
        "annotations_A2Inactive" = annotations_A2Inactive)
        )
}

groups <- getGroups(annotations = annotations)
annotations <- groups$annotations
annotations_A0A2 <- groups$annotations_A0A2
annotations_A1 <- groups$annotations_A1
annotations_A2 <- groups$annotations_A2
#annotations_A3 <- groups$annotations_A3
#annotations_A4 <- groups$annotations_A4
#annotations_AOld <- groups$annotations_AOld
#annotations_AOldER <- groups$annotations_AOldER
annotations_AVUS <- groups$annotations_AVUS
annotations_AVUSA2 <- groups$annotations_AVUSA2
annotations_A2Inactive <- groups$annotations_A2Inactive



# Function to recover probes in a given gene/Promoter
# @geneName:  name of the gene;
# @dataset: dataframe containing a column named "probe" with probes and a second columns with genes named "SYMBOL"
getProbesInGene <- function(geneName, dataset, which_column){
    if (missing(geneName)){
        stop("** Stop. Please provide a geneName to search for.")
    }
    if (missing(dataset)){
        stop("**Stop. Please provide a dataset to search in")
    }    
    if ("probe" %in% names(dataset) == FALSE){
        return("Please provide a dataframe with a column named 'probe'")
    }
    return(subset(dataset, dataset[[which_column]] == geneName)[, "probe"])
}



# Function returning a plot of annotations dataframes
# @annotations: annotation dataframe with samples in rows and containing all variables to plot as columns. ATTENTION, samples names must be provided as rownames
# @brewer: logical, use or not the RColorBrewer colors
# @forms: logical, add or not forms to tiles
PlotAnnotation <- function(annotations,
                           removedSamples  = NULL,
                           notUsedSamples  = NULL,
                           brewer          = TRUE,
                           sameColors      = FALSE,
                           forms           = FALSE,
                           transpose       = NULL)
    {
        require("pheatmap")
        require("ggplot2")
        require("ggsci")
        require("dplyr")
        require("tidyr")
        require("ggnewscale")
        # Modify the dataset
        annotations$rnames <- rownames(annotations)
        annotations %>%
            as_tibble() %>%
            mutate_all(as.factor) %>%
                   tidyr::gather(key = "key", value = "value", -rnames, factor_key = TRUE)  -> temp
        # Keep the order wanted
        Order.y <- as.vector(unique(temp$rnames))
        Order.x <- as.vector(rev(unique(temp$key)))
        # Relevel for the plot
        temp$rnames <- factor(temp$rnames, levels = Order.y)
        temp$key <- factor(temp$key, levels = Order.x)
        # Plot by adding a geaom_til for each variable
        plot <- ggplot()
        for (i in 1:length(unique(temp$key))){
            plot <- plot + 
                geom_tile(
                    data = temp %>%
                        filter(key==paste0(unique(temp$key)[i])) %>% droplevels,
                    aes(x = key, y = rnames, fill=value),
                    color = "white",
                    lwd = 1.5,
                    linetype = 1
                )
            # Add annotations is removed and not used samples were specified
            if (!is.null(removedSamples)){
                plot <- plot + geom_point(
                                   data = temp %>% filter(key==paste0(unique(temp$key)[i]) &
                                                          rnames %in% removedSamples) %>%
                                       droplevels,
                                   aes(x = key, y = rnames, shape = factor(1), color = factor(1)),
                                   size = 5, stroke = 3)
            }
            if (!is.null(notUsedSamples)){
                plot <- plot + geom_point(
                                   data = temp %>% filter(key==paste0(unique(temp$key)[i]) &
                                                          rnames %in% notUsedSamples) %>%
                                       droplevels,
                                   aes(x = key, y = rnames, shape = factor(2), color = factor(2)),
                                   size = 5, stroke = 3)
            }
            # Add the legend if removed and not used samples were specified
            if (!is.null(removedSamples) && !is.null(notUsedSamples)){
                plot <- plot + scale_shape_manual(name = "",
                                                  values = c(4, 4),
                                                  labels = c("Suppressed sample", "Not used sample")) +
                    scale_color_manual(name = "",
                                       values = c("red", "violet"),
                                       labels = c("Suppressed sample", "Not used sample"))
            }else{
                if (!is.null(removedSamples) && is.null(notUsedSamples)){
                    plot <- plot + scale_shape_manual(name = "",
                                                      values = 4,
                                                      labels = "Suppressed sample") +
                        scale_color_manual(name = "",
                                           values = c("red"),
                                           labels = c("Suppressed sample"))
                }
                if (is.null(removedSamples) && !is.null(notUsedSamples)){
                    plot <- plot + scale_shape_manual(name = "",
                                                      values = 4,
                                                      labels = "Not used sample") +
                        scale_color_manual(name = "",
                                           values = c("violet"),
                                           labels = c("Not used sample"))
                }
            }
            # Handle colors
            if (brewer == TRUE){
                require("RColorBrewer")
                pal <- data.frame(brewer.pal.info)
                pal <- subset(pal, pal$colorblind == "TRUE" & pal$category %in% c("div", "qual"))
                #print("step 1")
                pal <- rownames(pal)
                pal = rev(pal)
                #print("step 2")
                if (sameColors == FALSE){
                    #if (i > length(pal))
                    #{
                    #    if (i %% length(pal) == 0)
                    #    {
                    #        i = length(pal)
                    #    }else{
                    #        i = i %% 7
                    #    }
                    #}
                    #print(i)
                    plot <- plot + scale_fill_brewer(palette = pal[i],
                                                     guide = guide_legend(order = i),
                                                     name = paste0(unique(temp$key)[i]),
                                                     na.value = "grey")
                }else{
                    plot <- plot + scale_fill_brewer(palette = "Pastel2",
                                        #guide = guide_legend(order = i),
                                                     name = "",#paste0(unique(temp$key)[i]),
                                                     na.value = "grey")                
                }
            }else{
                plot <- plot + scale_fill_jco(name = paste0(unique(temp$key)[i]) )
            }
            plot <-  plot +
                scale_x_discrete(limits = Order.x) +
                    coord_fixed() +
                        ggnewscale::new_scale_fill()
        }
        plot <- plot + xlab("") + ylab("Samples") +
            theme_bw(base_size = 14) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  legend.position = "bottom",
                  legend.direction="vertical",
                  legend.margin=margin(),
                  legend.key.width = unit(0.7, "cm"),
                  legend.key.height = unit(0.7, "cm")
                  )
        if (!is.null(transpose))
        {
            plot <- plot + coord_flip()
        }
        return(plot)
    }



getM_C <- function(gMSet){
       Beta <- getBeta(gMSet)
       M <- log2(Beta / (1 - Beta))
       return(M)
}

getMCustom <- function(Beta){
       M <- log2(Beta / (1 - Beta))
       return(M)
}


variable_feature_selection <- function (data,
                                        thres.diff,
                                        thres.num,
                                        probs = 0.25
                                        )
{
    if (missing(thres.diff) && missing(thres.num))
    {
        stop("** Stop. No method found.")
    }
    if (!missing(thres.diff) && !missing(thres.num))
    {
        stop("** Stop. Choose one of the two options - thres.diff or thres.num")
    }
    diff.quantile <- function(x, probs = 1/length(x))
    {
        vv <- quantile(x, probs = c(probs, 1 - probs), na.rm = TRUE)
        return(vv[2] - vv[1])
    }
    rangeValues <- apply(data, 1, diff.quantile, probs = probs)
    if (!missing(thres.diff))
    {
        ind <- which(rangeValues >= thres.diff)
        genesList <- rownames(data[ind, ])
    }
    else if (!missing(thres.num))
    {
        genesList <- names(sort(rangeValues, decreasing = TRUE)[1:thres.num])
    }
    return(genesList)
}


create_dt <- function(x){
    DT::datatable(x,
                  filter = 'top',
                  extensions = 'Buttons',
                  options = list(dom = 'lfrtipB',
                                        # lengthMenu = 10,
                      lengthMenu = list(c(10,25,50,-1), c(10,25,50,"All")),
                      pageLength = 10,
                      autoWidth = FALSE,
                      pageLength=50,
                      scrollX='400px',
                      buttons = c('copy', 'excel'))
                  )
}

#dotplot_custom <- function(data, main=NULL){
#    to_plot <- data.frame(data@result[1:20, ])
#    to_plot <- fortify(to_plot, showCategory = 20, split=NULL)

#    to_plot <- dplyr::mutate(to_plot, x = eval(parse(text="GeneRatio")))
#    to_plot$GeneRatio2 <- to_plot$Count/sum(to_plot$Count)
#    idx <- order(to_plot$GeneRatio2, decreasing = TRUE)
#    to_plot$Description <- factor(to_plot$Description,
#                                  levels=rev(unique(to_plot$Description[idx])))
    
#    plot <- ggplot(to_plot, aes_string(x="GeneRatio2", y="Description", size="Count", color="p.adjust")) +
#        geom_point() + ggtitle(main) +
#            scale_color_continuous(low="red",
#                                   high="blue",
#                                   name = "p.adjust",
#                                   guide=guide_colorbar(reverse=TRUE)) +
#                                       scale_y_discrete(labels = as_labeller(30)) +
#                                           ylab(NULL) + xlab("GeneRatio") + DOSE::theme_dose(12) +
#                                               scale_size(range=c(3, 8)) + theme(plot.title = element_text(hjust = 0.5))
#    return(plot)
#}




# dotplot_custom is used to plot the enrichments of both GO and Kegg pathways. It takes as argument:
# @data: the dataframe returned by getEnrihedGO and getKeggEnriched
dotplot_custom <- function(data,
                           main              = NULL,
                           comparison        = FALSE,
                           orderAccordingTo  = "GeneRatio",
                           showCategory      = 45,
                           sign              = NULL)
{
    require(ggplot2)
    require(dplyr)
    require(stringr)
    require(readr)
    # Transform to dataframe if S4 object
    if (class(data) != "data.frame")
    {
        data <- data.frame(data)
    }

    if (dim(data)[1] < showCategory)
        {
            showCategory <- dim(data)[1]
        }
    if (orderAccordingTo  == "GeneRatio")
    {
        #to_plot <- dplyr::mutate(data, x = eval(parse(text=orderAccordingTo)))
        #print(to_plot)
        #to_plot$GeneRatio2 <- to_plot$Count/sum(to_plot$Count)
        to_plot <- data
        to_plot$Count <- to_plot$setSize
        to_plot$GeneRatio2 <- (stringr::str_count(to_plot$core_enrichment, "/") + 1) / to_plot$setSize
        #to_plot <- to_plot[order(to_plot$p.adjust, decreasing = FALSE), ]
        #to_plot <- to_plot[1:showCategory,]
        
        
        idx <- order(to_plot$GeneRatio2, decreasing = TRUE)
        to_plot <- to_plot[idx,]
        to_plot$Description <- factor(to_plot$Description,
                                      levels = rev(unique(to_plot$Description)))
        #to_plot$Description <- factor(to_plot$Description,
        #                              levels = rev(unique(to_plot$Description[idx])))
        #to_plot <- to_plot[idx,]
        to_plot <- to_plot[1:showCategory,]
        print(to_plot$Description)
        
        plot <- ggplot(to_plot, aes_string(x="GeneRatio2", y="Description", size="Count", color="p.adjust")) +
            geom_point() + ggtitle(main) +
            scale_color_continuous(low   = "red",
                                   high  = "blue",
                                   name  = "p.adjust",
                                   guide = guide_colorbar(reverse=TRUE)) +
            ##facet_grid(~comparison) +
            scale_y_discrete(labels = as_labeller(30)) +
            ylab(NULL) + xlab("GeneRatio") +
            DOSE::theme_dose(12) +
            scale_size(range=c(3, 8)) + theme(plot.title = element_text(hjust = 0.5))
    }else{
        if (orderAccordingTo == "signal")
        {
            data$signal <- stringr::str_split_fixed(data$leading_edge, pattern = "[,=]", n = 6)[,6]
            data$signal <- readr::parse_number(data$signal)/100
            orderAccordingTo <- "NES"
            size = "signal"
        }else{
            orderAccordingTo <- "NES"
            size = "NES"
        }
        to_plot <- dplyr::mutate(data, x = eval(parse(text=orderAccordingTo)))
        idx <- order(to_plot$NES, decreasing = TRUE)
        to_plot$Description <- factor(to_plot$Description,
                                         levels = rev(unique(to_plot$Description[idx])))
        #to_plot <- to_plot[abs(to_plot$NES) > 1.7,]
        plot <- ggplot(to_plot, aes_string(x = orderAccordingTo, y = "Description",
                                           size = size, color = "p.adjust")) +
            geom_point() + ggtitle(main) +
            scale_color_continuous(low   = "red",
                                   high  = "blue",
                                   name  = "p.adjust",
                                   guide = guide_colorbar(reverse=TRUE)) +
            ##facet_grid(~comparison) +
            scale_y_discrete(labels = as_labeller(30)) +
            ylab(NULL) + xlab("GeneRatio") +
            DOSE::theme_dose(12) +
            scale_size(range=c(3, 8)) + theme(plot.title = element_text(hjust = 0.5))
    }
    
    if (comparison != FALSE){
        plot <- plot +
            ggplot2::facet_grid(~comparison)
    }
    return(plot)
}







dotplot_custom_Paper<- function(data,
                           main              = NULL,
                           comparison        = FALSE,
                           orderAccordingTo  = "GeneRatio",
                           size              = NULL,
                           showSign          = NULL,
                           showCategory      = NULL,
                           sign              = NULL)
{
    require(ggplot2)
    require(dplyr)
    require(stringr)
    require(readr)
    # Transform to dataframe if S4 object
    if (class(data) != "data.frame")
    {
        data <- data.frame(data)
    }

    if (!is.null(showCategory)){
        if (dim(data)[1] < showCategory)
        {
            showCategory <- dim(data)[1]
        }
    }

    to_plot <- data
    to_plot$Count <- to_plot$setSize
    to_plot$GeneRatio2 <- (stringr::str_count(to_plot$core_enrichment, "/") + 1) / to_plot$setSize
        #to_plot <- to_plot[order(to_plot$p.adjust, decreasing = FALSE), ]
        #to_plot <- to_plot[1:showCategory,]
        
        
        idx <- order(to_plot$GeneRatio2, decreasing = TRUE)
        to_plot <- to_plot[idx,]
        to_plot$Description <- factor(to_plot$Description,
                                      levels = rev(unique(to_plot$Description)))
        #to_plot$Description <- factor(to_plot$Description,
        #                              levels = rev(unique(to_plot$Description[idx])))
        #to_plot <- to_plot[idx,]
        to_plot <- to_plot[1:showCategory,]
        print(to_plot$Description)



    

    
    to_plot <- data
    to_plot$Count <- to_plot$setSize

    if (orderAccordingTo == "GeneRatio"){
        to_plot$GeneRatio2 <- (stringr::str_count(to_plot$core_enrichment, "/") + 1) / to_plot$setSize
        orderAccordingTo = "GeneRatio2"
    }
#    idx <- order(to_plot[[orderAccordingTo]], decreasing = TRUE)
#    to_plot <- to_plot[idx,]

    if (showSign){
        #to_plot$Description <- ifelse(to_plot$NES < 0, paste0(to_plot$Description, " (-) "),
        #                              paste0(to_plot$Description, " (+) "))
        colText <- ifelse(to_plot$NES < 0, "#1C9E77", "#D95D00")
    }

    if (orderAccordingTo  == "NES" | size == "NES")
    {
        to_plot$NES <- abs(to_plot$NES)
        to_plot$Description <- factor(to_plot$Description,
                                      levels = rev(unique(to_plot$Description)))
    }
    
    if (orderAccordingTo == "signal" | size == "signal")
    {
        to_plot$signal <- stringr::str_split_fixed(to_plot$leading_edge, pattern = "[,=]", n = 6)[,6]
        to_plot$signal <- readr::parse_number(to_plot$signal)/100
    }
    
#    to_plot <- dplyr::mutate(to_plot, x = eval(parse(text=orderAccordingTo)))
    idx <- order(to_plot[[orderAccordingTo]], decreasing = TRUE)
    to_plot$Description <- factor(to_plot$Description,
                                  levels = rev(unique(to_plot$Description[idx])))

    # Keep only the n first categories if asked
    if (! is.null(showCategory))
    {
        to_plot <- to_plot[1:showCategory,]
        #print(to_plot$Description)
    }
    
    plot <- ggplot(to_plot, aes_string(x=orderAccordingTo, y="Description", size=size, color="p.adjust")) +
        geom_point() + ggtitle(main) +
        scale_color_continuous(low   = "red",
                               high  = "blue",
                               name  = "p.adjust",
                               guide = guide_colorbar(reverse=TRUE)) +
        scale_y_discrete(labels = as_labeller(30)) +
        ylab(NULL) + xlab("GeneRatio") +
        DOSE::theme_dose(12) +
        scale_size(range=c(3, 8)) + theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.y = element_text(colour = colText))
    
    if (comparison != FALSE){
        plot <- plot +
            ggplot2::facet_grid(~comparison)
    }
    return(plot)
}

#dotplot_custom(A1.Pro.KeggGse, orderAccordingTo = "GeneRatio", size = "NES", showSign=TRUE, showCategory = 50)








lolliplot_custom_Paper <- function(data,
                                   main              = NULL,
                                   comparison        = FALSE,
                                   orderAccordingTo  = "GeneRatio",
                                   size              = NULL,
                                   showSign          = NULL,
                                   showCategory      = NULL,
                                   sign              = NULL,
                                   categories        = NULL,
                                   ySize             = 10)
{
    require(ggplot2)
    require(dplyr)
    require(stringr)
    require(readr)

    print("here 1")
    #if (length(data) != 1){
    #    print("entered")
    #    if (class(data[[1]]) != "data.frame"){data[[1]] <- data.frame(data[[1]])}
    #    if (class(data[[2]]) != "data.frame"){data[[2]] <- data.frame(data[[2]])}
    #    data <- rbind(data[[1]], data[[2]])
    #} else{
    # Transform to dataframe if S4 object
        if (class(data) != "data.frame")
        {
            data <- data.frame(data)
        }
        
        if (!is.null(showCategory)){
            if (dim(data)[1] < showCategory)
            {
                showCategory <- dim(data)[1]
            }
        }
    #}
    print(data)
    to_plot <- data
    print(head(to_plot))
    to_plot$Count <- to_plot$setSize
    to_plot$GeneRatio2 <- (stringr::str_count(to_plot$core_enrichment, "/") + 1) / to_plot$setSize
        #to_plot <- to_plot[order(to_plot$p.adjust, decreasing = FALSE), ]
        #to_plot <- to_plot[1:showCategory,]

    to_plot$NES_sign <- ifelse(to_plot$NES <= 0, "hypomethylated", "hypermethylated")

    print(head(to_plot))
    #to_plot <- to_plot[1:showCategory,]
    
#    if (orderAccordingTo == "GeneRatio"){
#        to_plot$GeneRatio2 <- (stringr::str_count(to_plot$core_enrichment, "/") + 1) / to_plot$setSize
#        orderAccordingTo = "GeneRatio2"
#    }
#    idx <- order(to_plot[[orderAccordingTo]], decreasing = TRUE)
#    to_plot <- to_plot[idx,]

#    if (! is.null(showSign)){
        #to_plot$Description <- ifelse(to_plot$NES < 0, paste0(to_plot$Description, " (-) "),
        #                              paste0(to_plot$Description, " (+) "))
#        colText <- ifelse(to_plot$NES < 0, "#1C9E77", "#D95D00")
#    }

#    if (orderAccordingTo  == "NES" | size == "NES")
#    {
#        #to_plot$NES <- abs(to_plot$NES)
#        to_plot$Description <- factor(to_plot$Description,
#                                      levels = rev(unique(to_plot$Description)))
#    }
    
    if (orderAccordingTo == "signal" | size == "signal")
    {
        to_plot$signal <- stringr::str_split_fixed(to_plot$leading_edge, pattern = "[,=]", n = 6)[,6]
        to_plot$signal <- readr::parse_number(to_plot$signal)/100
    }
    
#    to_plot <- dplyr::mutate(to_plot, x = eval(parse(text=orderAccordingTo)))

    if (!comparison){
     #   idx <- order(to_plot$GeneRatio2, decreasing = TRUE)
     #   to_plot <- to_plot[idx,]
     #   to_plot$Description <- factor(to_plot$Description,
      #                                levels = rev(unique(to_plot$Description)))
            
            
        idx <- order(to_plot[["NES"]], decreasing = TRUE)
        to_plot$Description <- factor(to_plot$Description,
                                      levels = rev(unique(to_plot$Description[idx])))
    } else {
        #idx <- order(to_plot[to_plot$comparison == "A1", ]$NES, decreasing = TRUE)
        idx <- order(to_plot[["NES"]], decreasing = TRUE)
        to_plot$Description <- factor(to_plot$Description,
                                      levels = rev(unique(
                                          c(to_plot$Description[idx], to_plot[!to_plot$comparison == "A1", ]$Description)
                                      )))
    }
    # Keep only the n first categories if asked
    if (! is.null(showCategory))
    {
        to_plot <- to_plot[1:showCategory,]
    }

    #plot <- ggplot(to_plot, aes_string(x=orderAccordingTo, y="Description", size=size, color="p.adjust")) +
    #    geom_point() + ggtitle(main) +
    #    scale_color_continuous(low   = "red",
    #                           high  = "blue",
    #                           name  = "p.adjust",
    #                           guide = guide_colorbar(reverse=TRUE)) +
    #    scale_y_discrete(labels = as_labeller(30)) +
    #    ylab(NULL) + xlab("GeneRatio") +
    #    DOSE::theme_dose(12) +
    #    scale_size(range=c(3, 8)) + theme(plot.title = element_text(hjust = 0.5)) +
    #    theme(axis.text.y = element_text(colour = colText))

#    print(head(to_plot))
#    print(to_plot$NES)
#    print(to_plot[c("NES", "NES_sign")])
#    plot <- ggplot(to_plot, aes_stringx=orderAccordingTo, y="Description", size=size, color="p.adjust")) +

    med = mean(to_plot$p.adjust)
    med = min(to_plot$p.adjust) + ((max(to_plot$p.adjust) - min(to_plot$p.adjust))/2)
    require(scales)
    # This dataset `df_shorcut` and the aesthetics `aes(min, track)` are used
    # in the `geom_*` calls when they are not explicitly specified.
    #print(to_plot$Description)
    #if (! is.null(categories)){
    #    return(to_plot$Description[idx])
    #}
    plot <- ggplot(to_plot, aes(x=NES, y=Description, color = p.adjust)) +
          # Dotted line connection shortcut yes/no
          # This geom uses `df_connect` instead of `df_shorcut` because it is being
          # explicitly overridden
        geom_linerange(aes(xmin = 0, xmax = NES), color = "black", size = 1, linetype = "dotted") +
        geom_point(aes(x=NES, fill = p.adjust, size = GeneRatio2), colour = "black", pch=21) +
        geom_vline(xintercept = 0, linetype="dashed",
                   color = "black", size=1) +
        scale_y_discrete(labels = as_labeller(30)) +
        xlab("NES") + ylab("Pathway") + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = ySize)) +
        guides(size = guide_legend(override.aes=list(shape=21))) +
        labs(size = "GeneRatio") + 
        scale_fill_gradient2(low   = muted("red"),
                             high  = muted("blue"),
                             mid = "white",
                             name  = "adj. p-value",
                             midpoint = med, 
                             guide = guide_colorbar(reverse=TRUE))

    if (comparison){
        #label <- c("ATM PV vs all controls", "ATM PV ER+ vs ER+ controls", "A3", "A4", "ATM VUS vs all controls", "ATM VUS ER+ vs ER+ controls", "All ATM tumours vs all controls")
        #names(label) <- c("A1", "A2", "A3", "A4", "AVUS", "AVUSA2", "A0")
        label <- c("", "ATM PV ER+ vs ER+ controls", "ER+ ATM vs ER+ controls",
                   "ATM VUS ER+ vs ER+ controls", "ATM PV vs all controls",
                   "A3", "A4", "ATM VUS vs all controls", "All ATM tumours vs all controls", "ER+ ATM biallelic inactivation")
        names(label) <- c("comp", "A2", "A0A2",
                          "AVUSA2", "A1", "A3",
                          "A4", "AVUS", "A0", "A2Inactive")
        
        if (is.null(categories)){
            plot <- plot + #coord_flip() +
            #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            #facet_wrap(comparison ~ .)
                facet_wrap(. ~ comparison,
                           nrow = 1,
                           labeller = labeller(comparison = label)) #+
                   #scale_y_continuous(position = "right")

            print("I AM HERE")
        }

        if (! is.null(categories)){
            # Reorder labels
            Order <- ggplot_build(plot)$layout$panel_params[[1]]$y$get_labels()
            Order <- unlist(Order, recursive=FALSE)
            categories <- categories[Order, ]

            # Add the label annotations
            plot <- plot +
                ggnewscale::new_scale_fill() +
                geom_point(data = categories, aes(x=1,
                                                  y=c(1:dim(data.frame(categories))[1]),
                                                  colour = label, shape = label, group = label),
                           size=3) +
                scale_colour_manual(name = "categories", values=c("red","orange", "green", "blue","violet", "black")) +
                scale_shape_manual(name = "categories", values=c(3, 15, 19,17,18, 9))
            
            plot <- plot + facet_wrap(. ~ factor(comparison,
                                                 levels = c("comp", "A2", "A2Inactive", "A0A2",
                                                            "AVUSA2", "A1", "A3",
                                                            "A4", "AVUS", "A0")),
                                      nrow = 1,
                                      labeller = labeller(comparison = label))


            return(plot)
            require(grid)
            gt = ggplot_gtable(ggplot_build(plot))

            index <- gt$layout$l[grep('panel-1-1', gt$layout$name)]
            # modify the width
            gt$widths[index] = unit(0.01, "npc")#0.2*gt$widths[index]
            # Suppress x ticks
            index <- grep('axis-b-1-1', gt$layout$name)
            gt$grobs[[index]]$children[2]$axis[1]$grob[[1]]$gp$col <- "NA"
            gt$grobs[[index]]$children[2]$axis[1]$grob[[1]]$gp$fill <- "NA"
            gt$grobs[[index]]$children[2]$axis$grobs[[1]]$y <- rep(unit(0, units="cm"),5)
            gt$grobs[[index]]$children[2]$axis$grobs[[1]]$x <- rep(unit(0, units="cm"),5)

            index <- grep('panel-1-1', gt$layout$name)
            # Remove the line at 0
            val <- attributes(gt$grobs[[index]]$children[5])$names
            gt$grobs[[index]]$children[5][[val]]$gp$col <- "NA"

            # Remove the borders
            val <- attributes(gt$grobs[[index]]$children[8])$names
            gt$grobs[[index]]$children[8][[val]]$gp$col <- "NA"

            # Remove grid (x) #gTree[panel-1.gTree.24446]
            val <- attributes(gt$grobs[[index]]$children)$names[grep("grill.gTree", attributes(gt$grobs[[2]]$children)$names)]
            val2 <- attributes(gt$grobs[[index]]$children[[val]]$children)$names[grep("panel.grid.major.x", attributes(gt$grobs[[index]]$children[[val]]$children)$names)]
            gt$grobs[[index]]$children[[val]]$children[[val2]]$gp$col <- "NA"
            # Remove grid (y)
            val <- attributes(gt$grobs[[index]]$children)$names[grep("grill.gTree", attributes(gt$grobs[[index]]$children)$names)]
            val2 <- attributes(gt$grobs[[index]]$children[[val]]$children)$names[grep("panel.grid.major.y", attributes(gt$grobs[[index]]$children[[val]]$children)$names)]
            gt$grobs[[index]]$children[[val]]$children[[val2]]$gp$col <- "NA"

            # Remove header box of the categories # TableGrob strip
            index <- grep('strip-t-1-1', gt$layout$name)
            val <- attributes(gt$grobs[[index]]$grobs[[1]]$children)$names[1]
            gt$grobs[[index]]$grobs[[1]]$children[[val]]$gp$col <- "NA"
            gt$grobs[[index]]$grobs[[1]]$children[[val]]$gp$fill <- "NA"
            gt$grobs[[index]]$grobs[[1]]$children[[val]]$height <- unit(0, "npc")
            gt$grobs[[index]]$grobs[[1]]$children[[val]]$y <- unit(0, "npc")
            gt$grobs[[index]]$grobs[[1]]$children[[val]]$x <- unit(0, "npc")
            # remove text of the header
            val <- attributes(gt$grobs[[index]]$grobs[[1]]$children)$names[grep("text.x.top", attributes(gt$grobs[[index]]$grobs[[1]]$children)$names)]
            gt$grobs[[index]]$grobs[[1]]$children[[val]]$children[[1]]$label <- ""
            return(gt)
            
        }
    }
#+
#        DOSE::theme_dose(12)

                                        #scale_x_discrete(
#        geom_linerange(data = to_plot, 
#            aes(xmin = 0, xmax = NES, y = "Description"),
#            inherit.aes = FALSE,
#            color = "grey75",
#            linetype = "11" # dotted line
#        ) + 
                # Segment when shortcut==yes
                # When the `data` argument is missing in the `geom_*` function
#        geom_linerange(aes(xmin = min, xmax = NES, color = p.adjust), size = 2) +
        
                      # Segment when shortcut==no. Overlapped lineranges.
#        geom_linerange(data = df_longdist, aes(xmin = min, xmax = max, color = diff), size = 2) +
#        geom_linerange(data = df_longdist, aes(xmin = min, xmax = max), color = "#FFFCFC", size = .8)

    
    
 #   if (comparison != FALSE){
 #       plot <- plot +
 #           ggplot2::facet_grid(~comparison)
 #   }
    return(plot)
}










ggsaveDotPlot <- function(plot,
                          filename)
{
    Nrow <- nrow(ggplot_build(plot)$data[[1]])
    ggsave(plot = plot, filename = filename,
           height = Nrow * 0.25,
           width = 18,
           limitsize = FALSE)
}


VolcanoPlot <- function(dataset, pval_threshold, FC_threshold, comparison = ""){

    dat <- data.frame(FC = dataset$logFC, adj.P.value = dataset$adj.P.Val)
    dat$significant[dat$FC > FC_threshold & dat$adj.P.value < pval_threshold] <- "UP"
    dat$significant[dat$FC < - FC_threshold & dat$adj.P.value < pval_threshold] <- "DOWN"
    dat$significant[abs(dat$FC) < FC_threshold] <- "NOT SIGNIFICANT"
    dat$significant[abs(dat$FC) > FC_threshold & dat$adj.P.value > pval_threshold] <- "NOT SIGNIFICANT"
    
    NUp <- nrow(subset(dat, dat$significant == "UP"))
    NDown <- nrow(subset(dat, dat$significant == "DOWN"))
    
    #Visualization
    cols <- c("UP" = "red", "DOWN" = "blue", "NOT SIGNIFICANT" = "grey")
    VP1 <- ggplot(data=dat, aes(x=FC, y = -log10(adj.P.value), color=significant)) +
        geom_point(alpha=.6, size=1.2) +
        scale_colour_manual(values = cols) +
        geom_vline(xintercept = FC_threshold, colour="#990000", linetype="dashed") +
    geom_vline(xintercept = - FC_threshold, colour="#990000", linetype="dashed") +
    geom_hline(yintercept = -log10(pval_threshold), colour="#990000", linetype="dashed") +
    theme(legend.position="none") +
        xlab("log2(Fold Change)") +
        ylab("-log10 adjusted p-value") +
        theme_bw() +
        theme(legend.position = "none")
    VP_total <- ggdraw(add_sub(VP1, paste0(comparison, "\n", NUp + NDown, " probes differentially methylated (", NDown," downregulated and ", NUp, " upregulated)", "\nThresholds: ", FC_threshold, " FC; ", pval_threshold, " adj.p-value"), x = 0, hjust = 0, size = 11))
    
    return(VP_total)
}






assignGroups_func <- function(matrix,
                              sampleTable,
                              distance          = "pearson",
                              method            = "ward",
                              k                 = NULL,
                              #seed              = seed,
                              group_names       = paste0("group_",1:k),
                              group_num_column  = "group_num",
                              group_name_column = "group_name",
                              silhouette        = FALSE
                              )
{
    # Load packages
    require("ComplexHeatmap")
    require("reshape2")

    # Set seed
    #set.seed(seed)
    
    if (is.null(k))
    {
        stop("Please provide a number of group to assign (k)")
    }
    if (class(matrix) == "data.frame")
    {
        matrix = as.matrix(matrix)
    }
    # Make sure the length of group_names equals to k
    #if (!identical(x = integer(k), y = integer(length(group_names))))
    if (integer(k) != integer(length(group_names)))
    {
        stop("The length of group_names vector should be equal to k.",
             call.=TRUE
             )
    }
    # Hierarchical clustering
    pl <- Heatmap(mat                      = matrix,
                  clustering_distance_columns = distance,
                  clustering_method_columns   = method
                  )

    # Set seed
    #set.seed(seed)
    
    # Cut tree
    hc <- as.hclust(column_dend(pl))
    group_num <- cutree(tree = hc,
                        k    = k,
                        h    = NULL
                        )
    group <- melt(group_num)
    colnames(group) <- "group_num"
    
    # Add assigned group_name to sample table
    sampleTable$group_name <- group_names[group$group_num]

    # Silhouette
    if (silhouette)
    {
        # Compute silhouette
        silhouette <- silhouette_func(matrix,
                                      distance  = distance,
                                      method    = method,
                                      k         = k
                                      )
        # Add silhouette
        sampleTable$silhouette <- silhouette
    }
    return(sampleTable)
}



#################
# UMAP function #
#################


umap_func <- function(matrix,
                      sampleTable                = NULL,
                      title                      = "UMAP",
                      highlightedVar             = NULL,
                      highlightedVarName         = NULL,
                      colors                     = NULL,
                      n_neighbors                = 5,
                      metric                     = "euclidean",
                      seed                       = 9,
                      sampleLabel                = "Sample_Name",
                      sampleLabelColumnCondition = NULL,
                      sampleLabelValueCondition  = NULL,
                      fontsize                   = 10,
                      pointSize                  = 4,
                      stat_ellipse               = FALSE,
                      showLegend                 = TRUE,
                      x_axis_position            = "bottom",
                      assignGroup                = FALSE,
                      group_name_column          = c("ATM", "Non.ATM"),
                      k                          = 2
                      )
{

    # Load packages
    require(umap)
    require(ggplot2)
    require(ggrepel)
    
    # Set seed
    set.seed(seed)
    
    # Creating a projection
    umap <- umap(d = t(matrix),
                 n_neighbors = n_neighbors,
                 metric = metric
                 )

    # Extract main component
    umap_layout <- data.frame(umap$layout)
    colnames(umap_layout) <- c("UMAP1","UMAP2")
    
    # Assign groups if asked
    if (assignGroup)
    {
        sampleTable <- assignGroups_func(matrix,
                                         sampleTable       = sampleTable,
                                         distance          = metric,
                                         method            = "ward.D2",
                                         k                 = k,
                                         group_names       = paste0("group_",1:k),
                                         group_num_column  = "Sample_Group",
                                         group_name_column = group_name_column,
                                         silhouette        = FALSE
                                         )
        sampleTable$group_name <- factor(sampleTable$group_name)
    }

    # Combine PCA result with sample plan (if sample plan is provided)
    if (!is.null(sampleTable))
    {
        umap_layout_plus <- cbind(sampleTable,
                                  umap_layout
                                  )
        #print(head(umap_layout_plus))
    } else {
        umap_layout_plus <- umap_layout
    }

    print(head(umap_layout_plus))
    # Plot UMAP
    g <- ggplot(data = umap_layout_plus,
                mapping = aes(x = UMAP1, y = UMAP2)
                ) +
        scale_x_continuous(position = x_axis_position) +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        ggtitle(title) +
        theme_bw()+
        theme(text = element_text(size = fontsize),
              plot.title = element_text(hjust = 0.5))

    # Show legend or not
    if (showLegend)
    {
        g <- g + theme(legend.title = element_text(face = "bold"),
                       legend.key.height = unit(1.5,"cm")
                       )
    }
    else { g <- g + theme(legend.position = "none")}
    # Sample label
    if (!is.null(sampleLabel))
        {
            if (is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
                {
                    g <- g + geom_text_repel(mapping = aes(label = eval(parse(text=sampleLabel))),
                                             size    = fontsize/3)
                }
            else if (!is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
                {
                    stop("missing sampleLabelValueCondition, it is required since sampleLabelColumnCondition is provided", call.=TRUE)
                }
            else if (!is.null(sampleLabelColumnCondition) && !is.null(sampleLabelValueCondition))
                {
                    if (!is.element(el = sampleLabelValueCondition, set = as.vector(umap_layout_plus[,sampleLabelColumnCondition])))
                        {
                            warning(paste(sampleLabelColumnCondition,"column does not contain",sampleLabelValueCondition,"value. Therefore, no sample is labeled", sep = " "), call.=TRUE)
                        }
                    else {
                        g <- g +
                            geom_text_repel(data = dplyr::filter(.data = umap_layout_plus,
                                                eval(parse(text=sampleLabelColumnCondition)) == sampleLabelValueCondition
                                                                 ),
                                            mapping = aes(label = eval(parse(text=sampleLabel))),
                                            size = fontsize/5
                                            )
                    }
                }
        }

    # Highlight variable
    if (!is.null(highlightedVarName))
        {
        # Check if highlightedVar is a valid column of sample plan
        #if (!is.element(highlightedVar,sampleTable))
        #{
        #    stop("highlightedVar should be a valid column in sampleTable", call.=TRUE)
        #}
            if (assignGroup) {
                g <- g + geom_point(aes(color = eval(parse(text=highlightedVar)),
                                        shape = group_name),
                                    fill = "black",
                                    size = pointSize
                                    )
            }
            else {
                g <- g + geom_point(aes(color = eval(parse(text=highlightedVar))),
                                    fill = "black",
                                    size = pointSize
                                    ) + guides(color=guide_legend(title=highlightedVarName))
            }
        
            if (!is.null(colors)) {
                colors <- colors[1:length(unique(umap_layout_plus[[highlightedVar]]))]
                if (!is.null(names(colors))) {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors[levels(factor(umap_layout_plus[,highlightedVar]))]
                                               )
                    g <- g + scale_color_manual(name   = highlightedVarName,
                                            values = colors[levels(factor(umap_layout_plus[,highlightedVar]))]
                                                )
                
                } else {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors
                                               )
                    
                    g <- g + scale_color_manual(name   = highlightedVarName,
                                                values = colors
                                                )
                }
            }
            g <- g + scale_fill_discrete(name=highlightedVarName)
        } else {
            g <- g + geom_point(size = pointSize,
                                fill = "blue",
                                pch = 21
                                )
        }
    
    # Add ellipse
    if (stat_ellipse)
    {
        g <- g + stat_ellipse()
    }
    return(g)
}







umap_func_publi <- function(matrix,
                            returnData                 = FALSE,
                            sampleTable                = NULL,
                            title                      = "UMAP",
                            highlightedVar             = NULL,
                            highlightedVarName         = NULL,
                            colors                     = NULL,
                            n_neighbors                = 5,
                            metric                     = "euclidean",
                            seed                       = 9,
                            sampleLabel                = NULL, #"Sample_Name",
                            sampleLabelColumnCondition = NULL,
                            sampleLabelValueCondition  = NULL,
                            shape                      = NULL,
                            fontsize                   = 10,
                            pointSize                  = 4,
                            stat_ellipse               = FALSE,
                            showLegend                 = TRUE,
                            x_axis_position            = "bottom",
                            assignGroup                = FALSE,
                            group_name_column          = c("ATM", "Non.ATM"),
                            k                          = 2
                            )
{
    # Load packages
    require(umap)
    require(ggplot2)
    require(ggrepel)
    
    # Set seed
    set.seed(seed)
    
    # Creating a projection
    umap <- umap(d = t(matrix),
                 n_neighbors = n_neighbors,
                 metric = metric#,
                 #labels=scale(sampleTable$Diagnosis_age)
                 )
    #return(umap)
    # Extract main component
    umap_layout <- data.frame(umap$layout)
    colnames(umap_layout) <- c("UMAP1","UMAP2")
    
    # Assign groups if asked
    if (assignGroup)
    {
        sampleTable <- assignGroups_func(matrix,
                                         sampleTable       = sampleTable,
                                         distance          = metric,
                                         method            = "ward.D2",
                                         k                 = k,
                                         group_names       = paste0("group_",1:k),
                                         group_num_column  = "Sample_Group",
                                         group_name_column = group_name_column,
                                         silhouette        = FALSE
                                         )
        sampleTable$group_name <- factor(sampleTable$group_name)
    }

    # Combine PCA result with sample plan (if sample plan is provided)
    if (!is.null(sampleTable))
    {
        umap_layout_plus <- cbind(sampleTable,
                                  umap_layout
                                  )
        #print(head(umap_layout_plus))
    } else {
        umap_layout_plus <- umap_layout
    }

    if (returnData)
    {
        return(umap_layout_plus)
    }
    # Plot UMAP
    g <- ggplot(data = umap_layout_plus,
                mapping = aes(x = UMAP1, y = UMAP2)
                ) +
        scale_x_continuous(position = x_axis_position) +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        ggtitle(title) +
        theme_bw()+
        theme(text = element_text(size = fontsize),
              plot.title = element_text(hjust = 0.5))

    # Show legend or not
    if (showLegend)
    {
        g <- g + theme(legend.title = element_text(face = "bold"),
                       legend.key.height = unit(1.5,"cm")
                       )
    }
    else { g <- g + theme(legend.position = "none")}

    # Highlight variable
    if (!is.null(highlightedVarName))
        {
        # Check if highlightedVar is a valid column of sample plan
        #if (!is.element(highlightedVar,sampleTable))
        #{
        #    stop("highlightedVar should be a valid column in sampleTable", call.=TRUE)
        #}
            if (assignGroup) {
                g <- g + geom_point(aes(color = eval(parse(text=highlightedVar)),
                                        shape = group_name),
                                    fill = "black",
                                    size = pointSize
                                    )
            }
            else {
                g <- g + geom_point(aes(color = eval(parse(text=highlightedVar))),
                                    size = pointSize
                                    )
                if (!is.null(sampleLabelColumnCondition)){
                    g <- g + geom_point(data = dplyr::filter(.data = umap_layout_plus,
                                                             eval(parse(text=sampleLabelColumnCondition)) == sampleLabelValueCondition),
                                        aes(shape = eval(parse(text=shape)),
                                        #fill = eval(parse(text=highlightedVar))
                                            ),
                                        size = pointSize,
                                        color = "black",
                                        stroke = 1) +
                        scale_shape_manual(values=c(21,22,23,24)) +
                        guides(color=guide_legend(title=highlightedVarName),
                               shape= guide_legend(title="ATM status"))
                }
            }

            if (!is.null(colors)) {
                if (!is.null(names(colors))) {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors[levels(factor(umap_layout_plus[,highlightedVar]))]
                                               ) 
                } else {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors
                                               )
                    
                    g <- g + scale_color_manual(name   = highlightedVarName,
                                                values = colors
                                                )
                }
            }
            g <- g + scale_fill_discrete(name=highlightedVarName)
        } else {
            g <- g + geom_point(size = pointSize,
                                fill = "blue",
                                pch = 21
                                )
        }

    # Sample label
    if (!is.null(sampleLabel))
        {
            if (is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
                {
                    g <- g + geom_text_repel(mapping = aes(label = eval(parse(text=sampleLabel))),
                                             size    = fontsize/3)
                }
            else if (!is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
                {
                    stop("missing sampleLabelValueCondition, it is required since sampleLabelColumnCondition is provided", call.=TRUE)
                }
            else if (!is.null(sampleLabelColumnCondition) && !is.null(sampleLabelValueCondition))
                {
                    if (!is.element(el = sampleLabelValueCondition, set = as.vector(umap_layout_plus[,sampleLabelColumnCondition])))
                        {
                            warning(paste(sampleLabelColumnCondition,"column does not contain",sampleLabelValueCondition,"value. Therefore, no sample is labeled", sep = " "), call.=TRUE)
                        }
                    else {
                        g <- g +
                            geom_text_repel(data = dplyr::filter(.data = umap_layout_plus,
                                                eval(parse(text=sampleLabelColumnCondition)) == sampleLabelValueCondition
                                                                 ),
                                            mapping = aes(label = eval(parse(text=sampleLabel))),
                                            size = fontsize/4
                                            )
                    }
                }
        }
    
    # Add ellipse
    if (stat_ellipse)
    {
        g <- g + stat_ellipse()
    }
    return(g)
}





umap_func_publiv2 <- function(matrix,
                            returnData                 = FALSE,
                            sampleTable                = NULL,
                            title                      = "UMAP",
                            highlightedVar             = NULL,
                            highlightedVarName         = NULL,
                            colors                     = NULL,
                            n_neighbors                = 5,
                            metric                     = "euclidean",
                            seed                       = 9,
                            sampleLabel                = NULL, #"Sample_Name",
                            sampleLabelColumnCondition = NULL,
                            sampleLabelValueCondition  = NULL,
                            colCond = NULL,
                            valCond = NULL,
                            shape                      = NULL,
                            fontsize                   = 10,
                            pointSize                  = 4,
                            stat_ellipse               = FALSE,
                            showLegend                 = TRUE,
                            x_axis_position            = "bottom",
                            assignGroup                = FALSE,
                            group_name_column          = c("ATM", "Non.ATM"),
                            k                          = 2
                            )
{
    # Load packages
    require(umap)
    require(ggplot2)
    require(ggrepel)
    require('RColorBrewer')
    # Set seed
    set.seed(seed)
    
    # Creating a projection
    umap <- umap(d = t(matrix),
                 n_neighbors = n_neighbors,
                 metric = metric#,
                 #labels=scale(sampleTable$Diagnosis_age)
                 )
    #return(umap)
    # Extract main component
    umap_layout <- data.frame(umap$layout)
    colnames(umap_layout) <- c("UMAP1","UMAP2")

    if (highlightedVar == "ER_status")
    {
        annotations[["ER_status"]] <- as.character(annotations[["ER_status"]])       
        annotations[["ER_status2"]] <- annotations[["ER_status"]]
        annotations[is.na(annotations[["ER_status2"]]), ][["ER_status2"]] <- "NA"
        annotations[annotations[["ER_status2"]] == "1", ][["ER_status2"]] <- "Positive"
        annotations[annotations[["ER_status2"]] == "0", ][["ER_status2"]] <- "Negative"
        annotations[annotations[["ER_status2"]] == "NA", ][["ER_status2"]] <- NA
        print(annotations$ER_status2)
    }

    # Assign groups if asked
    if (assignGroup)
    {
        sampleTable <- assignGroups_func(matrix,
                                         sampleTable       = sampleTable,
                                         distance          = metric,
                                         method            = "ward.D2",
                                         k                 = k,
                                         group_names       = paste0("group_",1:k),
                                         group_num_column  = "Sample_Group",
                                         group_name_column = group_name_column,
                                         silhouette        = FALSE
                                         )
        sampleTable$group_name <- factor(sampleTable$group_name)
    }

    # Combine PCA result with sample plan (if sample plan is provided)
    if (!is.null(sampleTable))
    {
        umap_layout_plus <- cbind(sampleTable,
                                  umap_layout
                                  )
        #print(head(umap_layout_plus))
    } else {
        umap_layout_plus <- umap_layout
    }

    if (returnData)
    {
        return(umap_layout_plus)
    }
    # Plot UMAP
    g <- ggplot(data = umap_layout_plus,
                mapping = aes(x = UMAP1, y = UMAP2)
                ) +
        scale_x_continuous(position = x_axis_position) +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        ggtitle(title) +
        theme_bw()+
        theme(text = element_text(size = fontsize),
              plot.title = element_text(hjust = 0.5))

    # Show legend or not
    if (showLegend)
    {
        g <- g + theme(legend.title = element_text(face = "bold"),
                       legend.key.height = unit(1.5,"cm")
                       )
    }
    else { g <- g + theme(legend.position = "none")}

    print("all good")
    # Highlight variable
    if (!is.null(highlightedVarName))
        {
        # Check if highlightedVar is a valid column of sample plan
        #if (!is.element(highlightedVar,sampleTable))
        #{
        #    stop("highlightedVar should be a valid column in sampleTable", call.=TRUE)
        #}
            if (assignGroup) {
                g <- g + geom_point(aes(color = eval(parse(text=highlightedVar)),
                                        shape = group_name),
                                    fill = "black",
                                    size = pointSize
                                    ) +
                    guides(color=guide_legend(title=highlightedVarName),
                           shape= guide_legend(title="ATM status"))
            }
            else {
                g <- g + geom_point(aes(color = eval(parse(text=highlightedVar)),
                                        shape = eval(parse(text=shape))),
                                    size = pointSize
                                    )
                if (!is.null(sampleLabelColumnCondition)){
                    g <- g + geom_point(data = dplyr::filter(.data = umap_layout_plus,
                                                             eval(parse(text=sampleLabelColumnCondition)) == sampleLabelValueCondition),
                                        aes(shape = eval(parse(text=shape)),
                                        #fill = eval(parse(text=highlightedVar))
                                            ),
                                        size = pointSize,
                                        color = "black",
                                        stroke = 1) +
                        scale_shape_manual(values=c(21,22,23,24)) +
                        guides(color=guide_legend(title=highlightedVarName),
                               shape= guide_legend(title="ATM status"))
                }
            }
            print('all good 2')

            if (!is.null(colors)) {
                if (!is.null(names(colors))) {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors[levels(factor(umap_layout_plus[,highlightedVar]))]
                                               ) +
                    guides(color=guide_legend(title=highlightedVarName),
                           shape= guide_legend(title="ATM status"))
                } else {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors
                                               )
                    
                    g <- g + scale_color_manual(name   = highlightedVarName,
                                                values = colors
                                                ) +
                    guides(color=guide_legend(title=highlightedVarName),
                           shape= guide_legend(title="ATM status"))
                }
            }

                        print('all good 3')
            
            g <- g + scale_fill_discrete(name=highlightedVarName)
        } else {
            g <- g + geom_point(size = pointSize,
                                fill = "blue",
                                pch = 21
                                )
        }

    print('all good 4')
    if (!is.null(colCond)){
        g <- g + geom_point(data = dplyr::filter(.data = umap_layout_plus,
                                                   eval(parse(text=colCond)) == valCond),
                              aes(shape = eval(parse(text=shape))#,
                                  #fill = NA
                                  ),
                              shape = 1,
                              size = pointSize,
                             # color = "black",
                              stroke = 1)
       
      }

    print('return')
    return(g)
}








umap_func_publiv3 <- function(matrix,
                            returnData                 = FALSE,
                            sampleTable                = NULL,
                            title                      = "UMAP",
                            highlightedVar             = NULL,
                            highlightedVarName         = NULL,
                            colors                     = NULL,
                            n_neighbors                = 5,
                            metric                     = "euclidean",
                            seed                       = 9,
                            sampleLabel                = NULL, #"Sample_Name",
                            sampleLabelColumnCondition = NULL,
                            sampleLabelValueCondition  = NULL,
                            colCond = NULL,
                            valCond = NULL,
                            shape                      = NULL,
                            fontsize                   = 10,
                            pointSize                  = 4,
                            stat_ellipse               = FALSE,
                            showLegend                 = TRUE,
                            x_axis_position            = "bottom",
                            assignGroup                = FALSE,
                            group_name_column          = c("ATM", "Non.ATM"),
                            k                          = 2
                            )
{
    # Load packages
    require(umap)
    require(ggplot2)
    require(ggrepel)
    require('RColorBrewer')
    # Set seed
    set.seed(seed)
    
    # Creating a projection
    umap <- umap(d = t(matrix),
                 n_neighbors = n_neighbors,
                 metric = metric,
                 n_components = 3#,
                 #labels=scale(sampleTable$Diagnosis_age)
                 )
    #return(umap)
    # Extract main component
    umap_layout <- data.frame(umap$layout)
    colnames(umap_layout) <- c("UMAP1","UMAP2", "UMAP3")

    if (highlightedVar == "ER_status")
    {
        annotations[["ER_status"]] <- as.character(annotations[["ER_status"]])       
        annotations[["ER_status2"]] <- annotations[["ER_status"]]
        annotations[is.na(annotations[["ER_status2"]]), ][["ER_status2"]] <- "NA"
        annotations[annotations[["ER_status2"]] == "1", ][["ER_status2"]] <- "Positive"
        annotations[annotations[["ER_status2"]] == "0", ][["ER_status2"]] <- "Negative"
        annotations[annotations[["ER_status2"]] == "NA", ][["ER_status2"]] <- NA
        print(annotations$ER_status2)
    }

    # Assign groups if asked
    if (assignGroup)
    {
        sampleTable <- assignGroups_func(matrix,
                                         sampleTable       = sampleTable,
                                         distance          = metric,
                                         method            = "ward.D2",
                                         k                 = k,
                                         group_names       = paste0("group_",1:k),
                                         group_num_column  = "Sample_Group",
                                         group_name_column = group_name_column,
                                         silhouette        = FALSE
                                         )
        sampleTable$group_name <- factor(sampleTable$group_name)
    }

    # Combine PCA result with sample plan (if sample plan is provided)
    if (!is.null(sampleTable))
    {
        umap_layout_plus <- cbind(sampleTable,
                                  umap_layout
                                  )
        #print(head(umap_layout_plus))
    } else {
        umap_layout_plus <- umap_layout
    }

    if (returnData)
    {
        return(umap_layout_plus)
    }
    # Plot UMAP
    g <- ggplot(data = umap_layout_plus,
                mapping = aes(x = UMAP1, y = UMAP2)
                ) +
        scale_x_continuous(position = x_axis_position) +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        ggtitle(title) +
        theme_bw()+
        theme(text = element_text(size = fontsize),
              plot.title = element_text(hjust = 0.5))

    # Show legend or not
    if (showLegend)
    {
        g <- g + theme(legend.title = element_text(face = "bold"),
                       legend.key.height = unit(1.5,"cm")
                       )
    }
    else { g <- g + theme(legend.position = "none")}

    print("all good")
    # Highlight variable
    if (!is.null(highlightedVarName))
        {
        # Check if highlightedVar is a valid column of sample plan
        #if (!is.element(highlightedVar,sampleTable))
        #{
        #    stop("highlightedVar should be a valid column in sampleTable", call.=TRUE)
        #}
            if (assignGroup) {
                g <- g + geom_point(aes(color = eval(parse(text=highlightedVar)),
                                        shape = group_name),
                                    fill = "black",
                                    size = pointSize
                                    ) +
                    guides(color=guide_legend(title=highlightedVarName),
                           shape= guide_legend(title="ATM status"))
            }
            else {
                g <- g + geom_point(aes(color = eval(parse(text=highlightedVar)),
                                        shape = eval(parse(text=shape))),
                                    size = pointSize
                                    )
                if (!is.null(sampleLabelColumnCondition)){
                    g <- g + geom_point(data = dplyr::filter(.data = umap_layout_plus,
                                                             eval(parse(text=sampleLabelColumnCondition)) == sampleLabelValueCondition),
                                        aes(shape = eval(parse(text=shape)),
                                        #fill = eval(parse(text=highlightedVar))
                                            ),
                                        size = pointSize,
                                        color = "black",
                                        stroke = 1) +
                        scale_shape_manual(values=c(21,22,23,24)) +
                        guides(color=guide_legend(title=highlightedVarName),
                               shape= guide_legend(title="ATM status"))
                }
            }
            print('all good 2')

            if (!is.null(colors)) {
                if (!is.null(names(colors))) {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors[levels(factor(umap_layout_plus[,highlightedVar]))]
                                               ) +
                    guides(color=guide_legend(title=highlightedVarName),
                           shape= guide_legend(title="ATM status"))
                } else {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors
                                               )
                    
                    g <- g + scale_color_manual(name   = highlightedVarName,
                                                values = colors
                                                ) +
                    guides(color=guide_legend(title=highlightedVarName),
                           shape= guide_legend(title="ATM status"))
                }
            }

                        print('all good 3')
            
            g <- g + scale_fill_discrete(name=highlightedVarName)
        } else {
            g <- g + geom_point(size = pointSize,
                                fill = "blue",
                                pch = 21
                                )
        }

    print('all good 4')
    if (!is.null(colCond)){
        g <- g + geom_point(data = dplyr::filter(.data = umap_layout_plus,
                                                   eval(parse(text=colCond)) == valCond),
                              aes(shape = eval(parse(text=shape))#,
                                  #fill = NA
                                  ),
                              shape = 1,
                              size = pointSize,
                             # color = "black",
                              stroke = 1)
       
      }

    print('return')
    return(g)
}





DifferentialAnalysis <- function(data,
                                 annotations,
                                 colInterest,
                                 analysisLevel, #Gene or GenomeWide
                                 listProbes       = NULL,
                                 onlyCoding       = FALSE,
                                 weights          = NULL,
                                 contrastToUse    = "ATM - Non.ATM",
                                 comparison_text  = "Comparison of ATM and Non.ATM cancers",
                                 FC.threshold     = 1,
                                 p.val.threshold  = 0.05,
                                 path,
                                 height           = 7,
                                 aspect_ratio     = 2,
                                 GO               = FALSE,
                                 KEGG             = FALSE,
                                 keepAll          = FALSE
                                 )
{
    require('limma')
    require('ggplot2')
    require('cowplot')
    if (dim(data)[1] == 0)
    {
        stop("Provide a valid dataframe of M values")
    }
    if (is.null(annotations))
    {
        stop("Provide a valid dataframe of annotation")
    }
    #subset the data to have only samples present in the annotations
    data <- data[, rownames(annotations)]
    
    if (is.null(colInterest))
    {
        stop("Specify the columns of interest for the design matrix")
    }else
    {
        if (!colInterest %in% colnames(annotations))
        {
            stop("The column of interest is not in the annotations")
        }
    }
    if (analysisLevel == "Gene" && is.null(listProbes))
    {
        stop("Your analysis is at the gene level, provide a list of probes present un the gene")
    }else{
        if(!is.null(listProbes))
        {
            data <- subset(data, rownames(data) %in% listProbes)
        }
    }

    # If choosen, subset the data on coding genes, and translated pseudogenes
    if(onlyCoding)
    {
        require("biomaRt")
        if (!exists("mart"))
        {
            ensembl <-  useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", mirror="asia")
            mart <- useDataset(dataset="hsapiens_gene_ensembl", mart=ensembl)    
        }
        
        functions <- getBM(attributes= c("hgnc_symbol", "gene_biotype"),
                           values=rownames(data_Promoters),
                           mart= mart)
        functions <- functions[functions$gene_biotype %in% c("protein_coding", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"),]
        # get unique genes remaining
        funcGenes <- unique(functions$hgnc_symbol)
        # subset data
        data <- subset(data, rownames(data) %in% funcGenes)
    }
    
    # this is the factor of interest
    dataType <- factor(annotations[[colInterest]])
    # use the above to create a design matrix
    design <- model.matrix(~0+dataType, data=annotations)
    colnames(design) <- c(levels(dataType))    
    # fit the linear model
    if (is.null(weights)){
        fit <- lmFit(data, design)
    } else {
        fit <- lmFit(data, design, weights=weights)
    }
    # create a contrast matrix for specific comparisons
    contMatrix <- do.call(makeContrasts, list(contrast = contrastToUse, levels = design))
    contMatrix
    # fit the contrasts
    fit2 <- contrasts.fit(fit, contMatrix)
    fit2 <- eBayes(fit2)

    DMPs <- topTable(fit2, num=Inf, coef = 1, adjust = "BH", sort.by = "logFC")
    print('contrast fitted') 

    ##### Volcano plot	
    volcano <- VolcanoPlot(DMPs, p.val.threshold, FC.threshold, comparison_text)
    ggsave(paste0(path, "VolcanoPlot.png"), volcano, height=height, width=height*aspect_ratio)
    print('Volcano done')
    ##### List of probes
    if (keepAll)
    {
        DMPs$"DifferentiallyMethylated" <- ifelse(DMPs$adj.P.Val <= p.val.threshold & abs(DMPs$logFC) >= FC.threshold, "TRUE", "FALSE")

        write.csv(DMPs, paste(path, "DMPs.csv", sep=""))
        print('list of probes done')
        DMPs <- subset(DMPs,
                       DMPs$adj.P.Val <= p.val.threshold & abs(DMPs$logFC) >= FC.threshold)
        cat("Using a threshold of ", p.val.threshold, " on adjusted p-value and ", FC.threshold, " on FoldChange, ", nrow(DMPs), " probes were found differentially methylated (", nrow(DMPs[DMPs$logFC < 0 ,  ]), " downregulated and ", nrow(DMPs[DMPs$logFC >= 0 ,  ]), " upregulated).\n Enrichment analysis were then performed using those different probes."
)   
    }
    else
    {
        DMPs <- subset(DMPs,
                       DMPs$adj.P.Val <= p.val.threshold & abs(DMPs$logFC) >= FC.threshold)
        write.csv(DMPs, paste(path, "DMPs.csv", sep=""))
        print('list of probes done')
        cat("Using a threshold of ", p.val.threshold, " on adjusted p-value and ", FC.threshold, " on FoldChange, ", nrow(DMPs), " probes were found differentially methylated (", nrow(DMPs[DMPs$logFC < 0 ,  ]), " downregulated and ", nrow(DMPs[DMPs$logFC >= 0 ,  ]), " upregulated).\n
Enrichment analysis were then performed using those different probes."
)   
    }
    
}




    #print(heat_annotations$ClinVar)
    #heat_annotations$ClinVar[heat_annotations$ClinVar %in% c("Pathogenic", "Likely payhogenic", "Pathogenic/Likely pathogenic")] <- "PV"
    #heat_annotations$ClinVar[heat_annotations$ClinVar == ""] <- " "
    #heat_annotations$ClinVar <- factor(heat_annotations$ClinVar,
    #                                   levels = c("PV", " "))
    #print(heat_annotations$ClinVar)
    
    
    #heat_annotations$Study <- factor(heat_annotations$Study,
    #                                 levels = c("CoFAT", "GENESIS", "ABCFR", "MCCS"))






heatmap_func <- function(data,
                         annotations                = annotations,
                         otherAnnotations           = NULL,
                         AnnotationRB1              = NULL,
                         subsetAnnotation           = NULL,
#                         annotation_row            = annotation_row,
                         scale                      = "row",
                         color                      = NULL,
                         breaks                     = NA,
                         columns                    = NULL,
                         filename                   = NULL,
                         cluster_rows               = TRUE,
                         cluster_cols               = TRUE,
                         clustering_distance        = "euclidean",
                         cutree_rows                = 1,
                         cutree_cols                = 1,
                         show_rownames              = TRUE,
                         fontsize_row               = 12,
                         fontsize_col               = 12, 
                         treeheight_row             = 0,
                         treeheight_col             = 25,
                         width                      = 15*(5/4),
                         height                     = 15,
                         main                       = "",
                         seed                       = 9
                         )
{
    # Set seed
    set.seed(seed)
    
    # Load packages
    require("pheatmap")
    require("RColorBrewer")
    # Add double deletion (homozygote ATM variant)
    annotations$Homozygote  <- ifelse(annotations$Sample_Name %in% c("T0072", "T0249"), "Homozygous", "Heterozygous")
    # Add ATM LOH annotation
    ATM.LOH.Yes <- c("T0220", "T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
    ATM.LOH.No <- c("T0120", "T0015", "T0009", "T0249", "T0247")
    annotations$ATM_LOH  <- ifelse(annotations$Sample_Name %in% ATM.LOH.Yes, "LOH",
                                   ifelse(annotations$Sample_Name %in% ATM.LOH.No, "No LOH", " ")
                                   )
    
    # Create dataframe
    heat_annotations <- data.frame(row.names     = rownames(annotations),
                                   Study         = annotations$Study,
                                   Sample_Group  = annotations$Sample_Group,
                                   ER_status     = as.character(annotations$ER_status),
                                   PR_status     = as.character(annotations$PR_status),
                                   Subtype       = annotations$subtype,
                                   Variant_type  = annotations$Variant_type_aggregated,
                                   Grade         = annotations$Grade,
                                   TimeAfterExtraction = annotations$DegradationTime_Year_cat,
                                   ATM_LOH       = annotations$ATM_LOH,
                                   ClinVar       = annotations$ClinVar
                                   )
#    heat_annotations$Variant_type[heat_annotations$Variant_type == ""] <- NA
    heat_annotations[is.na(heat_annotations)] <- " "
    heat_annotations$Subtype <- factor(heat_annotations$Subtype,
                                       levels = c("Luminal A", "Luminal B", "Luminal A/B", "Luminal/HER2", "HER2", "TNBC", " "))
    heat_annotations[heat_annotations$PR_status == "Uncertain", ]$PR_status <- " "
    heat_annotations$Variant_type[heat_annotations$Variant_type == ""] <- " "
    heat_annotations$Variant_type <- factor(heat_annotations$Variant_type,
                                            levels = c("MV", "LoF", " "))
    heat_annotations$Grade[heat_annotations$Grade == NA] <- ""
    heat_annotations$Grade[heat_annotations$Grade == ""] <- " "
    heat_annotations$Grade <- factor(heat_annotations$Grade,
                                     levels = c("I", "II", "III", " "))
    heat_annotations$TimeAfterExtraction[heat_annotations$TimeAfterExtraction == ""] <-  " "
    heat_annotations$TimeAfterExtraction <- factor(heat_annotations$TimeAfterExtraction,
                                                   levels = c("1st quartile", "2nd quartile", "3rd quartile", "4th quartile", " "))

    print(heat_annotations$ClinVar)
    heat_annotations$ClinVar[heat_annotations$ClinVar %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")] <- "PV"
    heat_annotations$ClinVar[heat_annotations$ClinVar %in% c("Uncertain significance", "Likely benign", "Benign", "Benign/Likely benign")] <- "VSI"
    heat_annotations$ClinVar[heat_annotations$ClinVar == ""] <- " "
    heat_annotations$ClinVar[is.na(heat_annotations$ClinVar)] <- " "

    if ("VSI" %in% heat_annotations$ClinVar){
        heat_annotations$ClinVar <- factor(heat_annotations$ClinVar,
                                           levels = c("PV", "VSI", " "))
    } else {
        heat_annotations$ClinVar <- factor(heat_annotations$ClinVar,
                                           levels = c("PV", " "))
    }

    heat_annotations$Study <- factor(heat_annotations$Study,
                                     levels = c("CoFAT", "GENESIS", "ABCFR", "MCCS"))
    heat_anno_col <- list(
        Study = c(CoFAT = "#B4C300", GENESIS = "#70AD47", ABCFR = "#00DAE0", MCCS = "#e0a800"),    #GSE100850  
        Sample_Group  = c(ATM = "#B4C300", Non.ATM = "#C0A6FF"),
        ER_status     = c("0" = "#E199FF", "1" = "#FF9289", " " = "white"),
        PR_status     = c("0" = "#6ED102", "1" = "#E0B400", " " = "white"),
        Subtype       = c("Luminal A" = "#AA8F66", "Luminal B" = "#ED9B40", "Luminal A/B" = "#FFEEDB", "Luminal/HER2" = "#61C9A8", "HER2" = "#475657", "TNBC" = "#4A6D7C", " " = "white"),
        Variant_type  = c(LoF = "#CBA328", MV = "#EAD637", " " = "white"),
        Grade = c("I" = "#FEE08B", "II" = "#F46D43", "III" = "#A51626", " " = "white"),
#        TimeAfterExtraction = c("1st quartile" = "#333C83", "2nd quartile" = "#F24A72", "3rd quartile" = "#FDAF75", "4th quartile" = "#EAEA7F", " " = NA)
#        TimeAfterExtraction = c("1st quartile" = "#DD3497", "2nd quartile" = "#FA9FB5", "3rd quartile" = "#7FCDBB", "4th quartile" = "#1D91C0", " " = NA)
        TimeAfterExtraction = c("1st quartile" = "#C0A6FF", "2nd quartile" = "#00DAE0", "3rd quartile" = "#FF81F2", "4th quartile" = "#51BFFF", " " = "white"),
        ATM_LOH = c("LOH" = "#63B7FF", "No LOH" = "#000080", " " = "white"),

        ClinVar = if ("VSI" %in% heat_annotations$ClinVar){
            ClinVar = c("PV" = "#fc7303", "VSI" = "#fcba03", " " = "white")
        } else {
            print("VSI not present")
            print(heat_annotations$ClinVar)
            ClinVar = c("PV" = "#fc7303", " " = "white", "NA" = "white")
        }
        
    )
    print("heat_anno_col")
    print(heat_anno_col)
    print("printed")
    # replace colors
#    colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(90)
#    colors = colors[30:length(colors)]

    if (!is.null(subsetAnnotation))
    {
        heat_anno_col <- heat_anno_col[subsetAnnotation]
        heat_annotations <- heat_annotations[subsetAnnotation]
    }

    if (! is.null(otherAnnotations))
        {            
    # Add ATM LOH information and cytogenetic informations
#            ATM.LOH.Yes <- c("T0220", "T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
#            ATM.LOH.No <- c("T0120", "T0015", "T0009", "T0249", "T0247")
            
#            annotations$ATM_LOH  <- ifelse(annotations$Sample_Name %in% ATM.LOH.Yes, "LOH",
#                                           ifelse(annotations$Sample_Name %in% ATM.LOH.No, "No LOH", NA)
#                                           )
            # TUSC8
            L13q13.3.q14.11.Loss.LOH <- c("T0220", "T0076", "T0123", "T0248", "T0077R", "T0072", "T0078")
            L13q13.3.q14.11.No <- c("T0099", "T0249")
            
            annotations$TUSC8.L13q13.3.q14.11  <- ifelse(annotations$Sample_Name %in% L13q13.3.q14.11.Loss.LOH, "Loss/LOH",
                                                         ifelse(annotations$Sample_Name %in% L13q13.3.q14.11.No, "No alteration", NA)
                                                         )
            # RNU6
            L13q31.3.q32.1.Loss.LOH <- c("T0220", "T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
            L13q31.3.q32.1.No <- c("T0249")

            annotations$RNU6.L13q31.3.q32.1 <- ifelse(annotations$Sample_Name %in% L13q31.3.q32.1.Loss.LOH, "Loss/LOH",
                                                         ifelse(annotations$Sample_Name %in% L13q31.3.q32.1.No, "No alteration", NA)
                                                         )

            # COX5BP6
            L13q32.2.q32.3.Loss.LOH <- c("T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
            L13q32.2.q32.3.No <- c("T0220", "T0249")
            
            annotations$COX5BP6.L13q32.2.q32.3 <- ifelse(annotations$Sample_Name %in% L13q32.2.q32.3.Loss.LOH, "Loss/LOH",
                                                      ifelse(annotations$Sample_Name %in% L13q32.2.q32.3.No, "No alteration", NA)
                                                      )

            # NFATC3
            L16q22.1.Loss.LOH <- c("T0076", "T0249", "T0248", "T0077R", "T0078")
            L16q22.1.Loss <- c("T0099", "T0072")
            L16q22.1.No <- c("T0123")
            L16q22.1.Gain <- c("T0220")
            
            annotations$NFATC3.L16q22.1 <-
                ifelse(annotations$Sample_Name %in% L16q22.1.Loss.LOH, "Loss/LOH",
                       ifelse(annotations$Sample_Name %in% L16q22.1.Loss, "Loss",
                              ifelse(annotations$Sample_Name %in% L16q22.1.No, "No alteration",
                                     ifelse(annotations$Sample_Name %in% L16q22.1.Gain, "Gain", NA)
                                     )
                              )
                       )
            
            # APOL6
            L22q12.3.q13.31.Loss <- c("T0220")
            L22q12.3.q13.31.Loss.LOH <- c("T0076", "T0249", "T0077R", "T0072", "T0078")
            L22q12.3.q13.31.Loss.LOH.Loss <- c("T0099")
            L22q12.3.q13.31.No <- c("T0123", "T0248")
            
            annotations$APOL6.L22q12.3.q13.31 <-
                ifelse(annotations$Sample_Name %in% L22q12.3.q13.31.Loss, "Loss",
                       ifelse(annotations$Sample_Name %in% L22q12.3.q13.31.Loss.LOH, "Loss/LOH",
                              ifelse(annotations$Sample_Name %in% L22q12.3.q13.31.Loss.LOH.Loss, "Loss/LOH or loss",
                                     ifelse(annotations$Sample_Name %in% L22q12.3.q13.31.No, "No alteration", NA)
                                     )
                              )
                       )

            heat_annotations <-  merge(heat_annotations,
                  annotations[, c("TUSC8.L13q13.3.q14.11", "RNU6.L13q31.3.q32.1", "COX5BP6.L13q32.2.q32.3", "NFATC3.L16q22.1", "APOL6.L22q12.3.q13.31")], by = 0
                  )
            rownames(heat_annotations) <- heat_annotations$Row.names
            heat_annotations <- heat_annotations[, !names(heat_annotations) == "Row.names"]
            heat_annotations[is.na(heat_annotations)] <- " "
           
            dic <-  c("Loss" = "#C1FFC1", "Los/LOH" = "#44C97F", "Loss/LOH or Loss" = "#2D8C57", "Gain" = "#FC0502", "No alteration" = "#C1C1C1", " " = NA)
            
            heat_anno_col2 <- list(
                #ATM_LOH = c("LOH" = "#63B7FF", "No LOH" = "#000080", " " = NA),
                TUSC8.L13q13.3.q14.11 = c("Loss/LOH" = "#44C97F", "No alteration" = "#C1C1C1", " " = NA),
                RNU6.L13q31.3.q32.1 = c("Loss/LOH" = "#44C97F", "No alteration" = "#C1C1C1", " " = NA),
                COX5BP6.L13q32.2.q32.3 = c("Loss" = "#C1FFC1", "Loss/LOH" = "#44C97F", "Loss/LOH or Loss" = "#2D8C57", "Gain" = "#FC0502", "No alteration" = "#C1C1C1", " " = NA),
                NFATC3.L16q22.1 = c("Loss" = "#C1FFC1", "Loss/LOH" = "#44C97F", "Loss/LOH or Loss" = "#2D8C57", "Gain" = "#FC0502", "No alteration" = "#C1C1C1", " " = NA),
                APOL6.L22q12.3.q13.31 = c("Loss" = "#C1FFC1", "Loss/LOH" = "#44C97F", "Loss/LOH or loss" = "#2D8C57", "No alteration" = "#C1C1C1", " " = NA)
            )

            heat_anno_col <- c(heat_anno_col, heat_anno_col2)
        }

        if (! is.null(AnnotationRB1))
        {   
            # TUSC8
            L13q14.11.q21.2.Loss.LOH <- c("T0099", "T0220", "T0076", "T0123", "T0248", "T0077R", "T0072", "T0078")
            L13q14.11.q21.2.No <- c("T0249")
            
            annotations$RB1.L13q14.11.q21.2  <- ifelse(annotations$Sample_Name %in% L13q14.11.q21.2.Loss.LOH, "Loss/LOH",
                                                         ifelse(annotations$Sample_Name %in% L13q14.11.q21.2.No, "No alteration", NA)
                                                       )
            heat_annotations <-  merge(heat_annotations,
                  annotations[ c("RB1.L13q14.11.q21.2")], by = 0
                  )
            
            rownames(heat_annotations) <- heat_annotations$Row.names
            heat_annotations <- heat_annotations[, !names(heat_annotations) == "Row.names"]
            heat_annotations[is.na(heat_annotations)] <- " "
           
            dic <-  c("Loss" = "#C1FFC1", "Los/LOH" = "#44C97F", "Loss/LOH or Loss" = "#2D8C57", "Gain" = "#FC0502", "No alteration" = "#C1C1C1", " " = NA)
            
            heat_anno_col2 <- list(
                RB1.L13q14.11.q21.2 = c("Loss/LOH" = "#44C97F", "No alteration" = "#C1C1C1", " " = NA))

            heat_anno_col <- c(heat_anno_col, heat_anno_col2)
            
        }

    if (is.null(color))
    {
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
    }
    
    if (!is.null(filename))
    {
        pheatmap::pheatmap(data,
                           scale           = scale,
                           color           = color,
                           breaks          = breaks,
                           cluster_rows    = TRUE,
                           cluster_cols    = TRUE,
                           clustering_distance_rows  = clustering_distance,
                           clustering_distance_cols  = clustering_distance,
                           cutree_rows     = cutree_rows,
                           cutree_cols     = cutree_cols,
                           show_rownames   = TRUE,
                           fontsize_row    = fontsize_row,
                           fontsize_col    = fontsize_col,
 #                         annotation_row  = annotation_row,
                           annotation_col  = heat_annotations,
                           annotation_colors = heat_anno_col,
                           treeheight_row  = treeheight_row,
                           treeheight_col  = treeheight_col,
                           filename        = filename,
                           width           = width,
                           height          = height,
                           na_col          = "grey90",
                           main            = main,
                           silent = TRUE
                           )
    }else{
        #heat_annotations$Variant_type <- as.character(heat_annotations$Variant_type)
        heat_annotations[heat_annotations == " "] <- "NA"
        heat_annotations$Variant_type <- as.character(heat_annotations$Variant_type)
        heat_annotations[is.na(heat_annotations)] <- "NA"

        heat_annotations$ATM_LOH <- factor(heat_annotations$ATM_LOH,
                                           levels = c("LOH", "No LOH", "NA"))
        heat_annotations$Variant_type <- factor(heat_annotations$Variant_type,
                                           levels = c("LoF", "MV", "NA"))
        heat_annotations$ER_status <- factor(heat_annotations$ER_status,
                                           levels = c("0", "1", "NA"))

#        heat_annotations$Variant_type <- as.character(heat_annotations$Variant_type)
#        heat_annotations[is.na(heat_annotations)] <- " "

 #       heat_annotations$ATM_LOH <- factor(heat_annotations$ATM_LOH,
 #                                          levels = c("LOH", "No LOH", " "))
 #       heat_annotations$Variant_type <- factor(heat_annotations$Variant_type,
 #                                          levels = c("LoF", "MV", " "))
        
        for (i in seq_along(names(heat_anno_col)))
        {
            print(names(heat_anno_col)[i])
            names(heat_anno_col[[names(heat_anno_col[i])]])[names(heat_anno_col[[names(heat_anno_col[i])]]) == " "] <- "NA"
        }

        print("Variant type")
        print(heat_annotations$Variant_type)
        print(heat_anno_col$Variant_type)
        
        print("ER status")
        print(heat_annotations$ER_status)
        print(heat_anno_col$ER_status)

        print("clinvar")
        print(heat_annotations$ClinVar)
        print(heat_anno_col$ClinVar)

        plot <- ComplexHeatmap::pheatmap(data,
                                         name = " ",
                                         row_split = NULL,
                                         scale            = scale,
                                         color            = color,
                                        #breaks           = breaks,
                                        #cluster_rows     = TRUE,
                                        #cluster_cols     = TRUE,
                                         clustering_distance_rows  = clustering_distance,
                                         clustering_distance_cols  = clustering_distance,
                                        #cutree_rows      = cutree_rows,
                                        #cutree_cols      = cutree_cols,
                                         show_rownames    = TRUE,
                                         fontsize_row     = fontsize_row,
                                         fontsize_col     = fontsize_col,
                                        #                          annotation_row  = annotation_row,
                                         annotation_col   = heat_annotations,
                                         annotation_colors = heat_anno_col,
                                         treeheight_row   = treeheight_row,
                                         treeheight_col   = treeheight_col,
                                        #width            = width,
                                        #height           = height,
                                         #na_col           = "grey90",
                                         na_col           = "#FFFFFF",
                                         main             = main,
                                        #silent           = TRUE
                                         heatmap_legend_param = list(legend_height = unit(4, "cm")),
                                         show_row_dend  = FALSE
                                       #, border = "#98999C")
                                         )
        #print(ComplexHeatmap::list_components())
#        print(grid::seekViewport())
#        print(ComplexHeatmap::seekViewport())

        co = ComplexHeatmap::column_order(
                                 ComplexHeatmap::draw(plot,
                                                      annotation_legend_list =
                                                          ComplexHeatmap::Legend(pch = 17, legend_gp = grid::gpar(col = "#AE217E"), type = "points", labels = "Homozygous")))
        pos <- which(annotations[co, c("Sample_Name", "Homozygote")]$Homozygote == "Homozygous")
        plot <- plot + ComplexHeatmap::decorate_annotation("ATM_LOH",
        {
            od = co[[1]]
            print(od)
            grid::grid.points(x = pos,
                              y = rep(unit(c(0.5), "npc"), length(pos)),
                              pch = 17,
                              gp = grid::gpar(col = "#AE217E")
                              )
        })
        return(plot)
        #return(ComplexHeatmap::draw(plot, merge_legends = TRUE, legend_grouping = "original"))
    }
}






























########################################################################
################# Duplicate unmodified #################################
########################################################################
#
# Function to plot a heatmap with annotations
# @data:
# @annotations

heatmap_func_first<- function(data,
                         annotations                = annotations,
                         otherAnnotations           = NULL,
                         AnnotationRB1              = NULL,
                         subsetAnnotation           = NULL,
#                         annotation_row            = annotation_row,
                         scale                      = "row",
                         color                      = NULL,
                         breaks                     = NA,
                         columns                    = NULL,
                         filename                   = NULL,
                         cluster_rows               = TRUE,
                         cluster_cols               = TRUE,
                         clustering_distance        = "euclidean",
                         cutree_rows                = 1,
                         cutree_cols                = 1,
                         show_rownames              = TRUE,
                         fontsize_row               = 12,
                         fontsize_col               = 12, 
                         treeheight_row             = 0,
                         treeheight_col             = 25,
                         width                      = 15*(5/4),
                         height                     = 15,
                         main                       = "",
                         seed                       = 9
                         )
{
    # Set seed
    set.seed(seed)
    
    # Load packages
    require("pheatmap")
    require("RColorBrewer")
    # Add double deletion (homozygote ATM variant)
    annotations$Homozygote  <- ifelse(annotations$Sample_Name %in% c("T0072", "T0249"), "Homozygous", "Heterozygous")
    # Add ATM LOH annotation
    ATM.LOH.Yes <- c("T0220", "T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
    ATM.LOH.No <- c("T0120", "T0015", "T0009", "T0249", "T0247")
    annotations$ATM_LOH  <- ifelse(annotations$Sample_Name %in% ATM.LOH.Yes, "LOH",
                                   ifelse(annotations$Sample_Name %in% ATM.LOH.No, "No LOH", " ")
                                   )
    
    # Create dataframe
    heat_annotations <- data.frame(row.names     = rownames(annotations),
                                   Study         = annotations$Study,
                                   Sample_Group  = annotations$Sample_Group,
                                   ER_status     = as.character(annotations$ER_status),
                                   PR_status     = as.character(annotations$PR_status),
                                   Subtype       = annotations$subtype,
                                   Variant_type  = annotations$Variant_type_aggregated,
                                   Grade         = annotations$Grade,
                                   TimeAfterExtraction = annotations$DegradationTime_Year_cat,
                                   ATM_LOH       = annotations$ATM_LOH
                                   )
#    heat_annotations$Variant_type[heat_annotations$Variant_type == ""] <- NA
    heat_annotations[is.na(heat_annotations)] <- " "
    heat_annotations$Subtype <- factor(heat_annotations$Subtype,
                                       levels = c("Luminal A", "Luminal B", "Luminal A/B", "Luminal/HER2", "HER2", "TNBC", " "))
    heat_annotations$Variant_type[heat_annotations$Variant_type == ""] <- " "
    heat_annotations$Variant_type <- factor(heat_annotations$Variant_type,
                                            levels = c("MV", "LoF", " "))
    heat_annotations$Grade[heat_annotations$Grade == NA] <- ""
    heat_annotations$Grade[heat_annotations$Grade == ""] <- " "
    heat_annotations$Grade <- factor(heat_annotations$Grade,
                                     levels = c("I", "II", "III", " "))
    heat_annotations$TimeAfterExtraction[heat_annotations$TimeAfterExtraction == ""] <-  " "
    heat_annotations$TimeAfterExtraction <- factor(heat_annotations$TimeAfterExtraction,
                                                   levels = c("1st quartile", "2nd quartile", "3rd quartile", "4th quartile", " "))
    heat_anno_col <- list(
        Study = c(CoFAT = "#B4C300", GENESIS = "#70AD47", ABCFR = "#00DAE0", MCCS = "#e0a800"),    #GSE100850  
        Sample_Group  = c(ATM = "#B4C300", Non.ATM = "#C0A6FF"),
        ER_status     = c("0" = "#E199FF", "1" = "#FF9289", " " = "white"),
        PR_status     = c("0" = "#6ED102", "1" = "#E0B400", " " = "white"),
        Subtype       = c("Luminal A" = "#AA8F66", "Luminal B" = "#ED9B40", "Luminal A/B" = "#FFEEDB", "Luminal/HER2" = "#61C9A8", "HER2" = "#475657", "TNBC" = "#4A6D7C", " " = "white"),
        Variant_type  = c(LoF = "#CBA328", MV = "#EAD637", " " = "white"),
        Grade = c("I" = "#FEE08B", "II" = "#F46D43", "III" = "#A51626", " " = "white"),
#        TimeAfterExtraction = c("1st quartile" = "#333C83", "2nd quartile" = "#F24A72", "3rd quartile" = "#FDAF75", "4th quartile" = "#EAEA7F", " " = NA)
#        TimeAfterExtraction = c("1st quartile" = "#DD3497", "2nd quartile" = "#FA9FB5", "3rd quartile" = "#7FCDBB", "4th quartile" = "#1D91C0", " " = NA)
        TimeAfterExtraction = c("1st quartile" = "#C0A6FF", "2nd quartile" = "#00DAE0", "3rd quartile" = "#FF81F2", "4th quartile" = "#51BFFF", " " = "white"),
        ATM_LOH = c("LOH" = "#63B7FF", "No LOH" = "#000080", " " = "white")
    )
    # replace colors
#    colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(90)
#    colors = colors[30:length(colors)]

    if (!is.null(subsetAnnotation))
    {
        heat_anno_col <- heat_anno_col[subsetAnnotation]
        heat_annotations <- heat_annotations[subsetAnnotation]
    }

    if (! is.null(otherAnnotations))
        {            
    # Add ATM LOH information and cytogenetic informations
#            ATM.LOH.Yes <- c("T0220", "T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
#            ATM.LOH.No <- c("T0120", "T0015", "T0009", "T0249", "T0247")
            
#            annotations$ATM_LOH  <- ifelse(annotations$Sample_Name %in% ATM.LOH.Yes, "LOH",
#                                           ifelse(annotations$Sample_Name %in% ATM.LOH.No, "No LOH", NA)
#                                           )
            # TUSC8
            L13q13.3.q14.11.Loss.LOH <- c("T0220", "T0076", "T0123", "T0248", "T0077R", "T0072", "T0078")
            L13q13.3.q14.11.No <- c("T0099", "T0249")
            
            annotations$TUSC8.L13q13.3.q14.11  <- ifelse(annotations$Sample_Name %in% L13q13.3.q14.11.Loss.LOH, "Loss/LOH",
                                                         ifelse(annotations$Sample_Name %in% L13q13.3.q14.11.No, "No alteration", NA)
                                                         )
            # RNU6
            L13q31.3.q32.1.Loss.LOH <- c("T0220", "T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
            L13q31.3.q32.1.No <- c("T0249")

            annotations$RNU6.L13q31.3.q32.1 <- ifelse(annotations$Sample_Name %in% L13q31.3.q32.1.Loss.LOH, "Loss/LOH",
                                                         ifelse(annotations$Sample_Name %in% L13q31.3.q32.1.No, "No alteration", NA)
                                                         )

            # COX5BP6
            L13q32.2.q32.3.Loss.LOH <- c("T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
            L13q32.2.q32.3.No <- c("T0220", "T0249")
            
            annotations$COX5BP6.L13q32.2.q32.3 <- ifelse(annotations$Sample_Name %in% L13q32.2.q32.3.Loss.LOH, "Loss/LOH",
                                                      ifelse(annotations$Sample_Name %in% L13q32.2.q32.3.No, "No alteration", NA)
                                                      )

            # NFATC3
            L16q22.1.Loss.LOH <- c("T0076", "T0249", "T0248", "T0077R", "T0078")
            L16q22.1.Loss <- c("T0099", "T0072")
            L16q22.1.No <- c("T0123")
            L16q22.1.Gain <- c("T0220")
            
            annotations$NFATC3.L16q22.1 <-
                ifelse(annotations$Sample_Name %in% L16q22.1.Loss.LOH, "Loss/LOH",
                       ifelse(annotations$Sample_Name %in% L16q22.1.Loss, "Loss",
                              ifelse(annotations$Sample_Name %in% L16q22.1.No, "No alteration",
                                     ifelse(annotations$Sample_Name %in% L16q22.1.Gain, "Gain", NA)
                                     )
                              )
                       )
            
            # APOL6
            L22q12.3.q13.31.Loss <- c("T0220")
            L22q12.3.q13.31.Loss.LOH <- c("T0076", "T0249", "T0077R", "T0072", "T0078")
            L22q12.3.q13.31.Loss.LOH.Loss <- c("T0099")
            L22q12.3.q13.31.No <- c("T0123", "T0248")
            
            annotations$APOL6.L22q12.3.q13.31 <-
                ifelse(annotations$Sample_Name %in% L22q12.3.q13.31.Loss, "Loss",
                       ifelse(annotations$Sample_Name %in% L22q12.3.q13.31.Loss.LOH, "Loss/LOH",
                              ifelse(annotations$Sample_Name %in% L22q12.3.q13.31.Loss.LOH.Loss, "Loss/LOH or loss",
                                     ifelse(annotations$Sample_Name %in% L22q12.3.q13.31.No, "No alteration", NA)
                                     )
                              )
                       )

            heat_annotations <-  merge(heat_annotations,
                  annotations[, c("TUSC8.L13q13.3.q14.11", "RNU6.L13q31.3.q32.1", "COX5BP6.L13q32.2.q32.3", "NFATC3.L16q22.1", "APOL6.L22q12.3.q13.31")], by = 0
                  )
            rownames(heat_annotations) <- heat_annotations$Row.names
            heat_annotations <- heat_annotations[, !names(heat_annotations) == "Row.names"]
            heat_annotations[is.na(heat_annotations)] <- " "
           
            dic <-  c("Loss" = "#C1FFC1", "Los/LOH" = "#44C97F", "Loss/LOH or Loss" = "#2D8C57", "Gain" = "#FC0502", "No alteration" = "#C1C1C1", " " = NA)
            
            heat_anno_col2 <- list(
                #ATM_LOH = c("LOH" = "#63B7FF", "No LOH" = "#000080", " " = NA),
                TUSC8.L13q13.3.q14.11 = c("Loss/LOH" = "#44C97F", "No alteration" = "#C1C1C1", " " = NA),
                RNU6.L13q31.3.q32.1 = c("Loss/LOH" = "#44C97F", "No alteration" = "#C1C1C1", " " = NA),
                COX5BP6.L13q32.2.q32.3 = c("Loss" = "#C1FFC1", "Loss/LOH" = "#44C97F", "Loss/LOH or Loss" = "#2D8C57", "Gain" = "#FC0502", "No alteration" = "#C1C1C1", " " = NA),
                NFATC3.L16q22.1 = c("Loss" = "#C1FFC1", "Loss/LOH" = "#44C97F", "Loss/LOH or Loss" = "#2D8C57", "Gain" = "#FC0502", "No alteration" = "#C1C1C1", " " = NA),
                APOL6.L22q12.3.q13.31 = c("Loss" = "#C1FFC1", "Loss/LOH" = "#44C97F", "Loss/LOH or loss" = "#2D8C57", "No alteration" = "#C1C1C1", " " = NA)
            )

            heat_anno_col <- c(heat_anno_col, heat_anno_col2)
        }

        if (! is.null(AnnotationRB1))
        {   
            # TUSC8
            L13q14.11.q21.2.Loss.LOH <- c("T0099", "T0220", "T0076", "T0123", "T0248", "T0077R", "T0072", "T0078")
            L13q14.11.q21.2.No <- c("T0249")
            
            annotations$RB1.L13q14.11.q21.2  <- ifelse(annotations$Sample_Name %in% L13q14.11.q21.2.Loss.LOH, "Loss/LOH",
                                                         ifelse(annotations$Sample_Name %in% L13q14.11.q21.2.No, "No alteration", NA)
                                                       )
            heat_annotations <-  merge(heat_annotations,
                  annotations[ c("RB1.L13q14.11.q21.2")], by = 0
                  )
            
            rownames(heat_annotations) <- heat_annotations$Row.names
            heat_annotations <- heat_annotations[, !names(heat_annotations) == "Row.names"]
            heat_annotations[is.na(heat_annotations)] <- " "
           
            dic <-  c("Loss" = "#C1FFC1", "Los/LOH" = "#44C97F", "Loss/LOH or Loss" = "#2D8C57", "Gain" = "#FC0502", "No alteration" = "#C1C1C1", " " = NA)
            
            heat_anno_col2 <- list(
                RB1.L13q14.11.q21.2 = c("Loss/LOH" = "#44C97F", "No alteration" = "#C1C1C1", " " = NA))

            heat_anno_col <- c(heat_anno_col, heat_anno_col2)
            
        }

    if (is.null(color))
    {
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
    }
    
    if (!is.null(filename))
    {
        pheatmap::pheatmap(data,
                           scale           = scale,
                           color           = color,
                           breaks          = breaks,
                           cluster_rows    = TRUE,
                           cluster_cols    = TRUE,
                           clustering_distance_rows  = clustering_distance,
                           clustering_distance_cols  = clustering_distance,
                           cutree_rows     = cutree_rows,
                           cutree_cols     = cutree_cols,
                           show_rownames   = TRUE,
                           fontsize_row    = fontsize_row,
                           fontsize_col    = fontsize_col,
 #                         annotation_row  = annotation_row,
                           annotation_col  = heat_annotations,
                           annotation_colors = heat_anno_col,
                           treeheight_row  = treeheight_row,
                           treeheight_col  = treeheight_col,
                           filename        = filename,
                           width           = width,
                           height          = height,
                           na_col          = "grey90",
                           main            = main,
                           silent = TRUE
                           )
    }else{
        #heat_annotations$Variant_type <- as.character(heat_annotations$Variant_type)
        heat_annotations[heat_annotations == " "] <- "NA"
        heat_annotations$Variant_type <- as.character(heat_annotations$Variant_type)
        heat_annotations[is.na(heat_annotations)] <- "NA"

        heat_annotations$ATM_LOH <- factor(heat_annotations$ATM_LOH,
                                           levels = c("LOH", "No LOH", "NA"))
        heat_annotations$Variant_type <- factor(heat_annotations$Variant_type,
                                           levels = c("LoF", "MV", "NA"))
        heat_annotations$ER_status <- factor(heat_annotations$ER_status,
                                           levels = c("0", "1", "NA"))

#        heat_annotations$Variant_type <- as.character(heat_annotations$Variant_type)
#        heat_annotations[is.na(heat_annotations)] <- " "

 #       heat_annotations$ATM_LOH <- factor(heat_annotations$ATM_LOH,
 #                                          levels = c("LOH", "No LOH", " "))
 #       heat_annotations$Variant_type <- factor(heat_annotations$Variant_type,
 #                                          levels = c("LoF", "MV", " "))
        
        for (i in seq_along(names(heat_anno_col)))
        {
            names(heat_anno_col[[names(heat_anno_col[i])]])[names(heat_anno_col[[names(heat_anno_col[i])]]) == " "] <- "NA"
        }

        print("Variant type")
        print(heat_annotations$Variant_type)
        print(heat_anno_col$Variant_type)
        
        print("ER status")
        print(heat_annotations$ER_status)
        print(heat_anno_col$ER_status)

        plot <- ComplexHeatmap::pheatmap(data,
                                         name = " ",
                                         row_split = NULL,
                                         scale            = scale,
                                         color            = color,
                                        #breaks           = breaks,
                                        #cluster_rows     = TRUE,
                                        #cluster_cols     = TRUE,
                                         clustering_distance_rows  = clustering_distance,
                                         clustering_distance_cols  = clustering_distance,
                                        #cutree_rows      = cutree_rows,
                                        #cutree_cols      = cutree_cols,
                                         show_rownames    = TRUE,
                                         fontsize_row     = fontsize_row,
                                         fontsize_col     = fontsize_col,
                                        #                          annotation_row  = annotation_row,
                                         annotation_col   = heat_annotations,
                                         annotation_colors = heat_anno_col,
                                         treeheight_row   = treeheight_row,
                                         treeheight_col   = treeheight_col,
                                        #width            = width,
                                        #height           = height,
                                         #na_col           = "grey90",
                                         na_col           = "#FFFFFF",
                                         main             = main,
                                        #silent           = TRUE
                                         heatmap_legend_param = list(legend_height = unit(4, "cm")),
                                         show_row_dend  = FALSE
                                       #, border = "#98999C")
                                         )
        #print(ComplexHeatmap::list_components())
#        print(grid::seekViewport())
#        print(ComplexHeatmap::seekViewport())

        co = ComplexHeatmap::column_order(
                                 ComplexHeatmap::draw(plot,
                                                      annotation_legend_list =
                                                          ComplexHeatmap::Legend(pch = 17, legend_gp = grid::gpar(col = "#AE217E"), type = "points", labels = "Homozygous")))
        pos <- which(annotations[co, c("Sample_Name", "Homozygote")]$Homozygote == "Homozygous")
        plot <- plot + ComplexHeatmap::decorate_annotation("ATM_LOH",
        {
            od = co[[1]]
            print(od)
            grid::grid.points(x = pos,
                              y = rep(unit(c(0.5), "npc"), length(pos)),
                              pch = 17,
                              gp = grid::gpar(col = "#AE217E")
                              )
        })
        return(plot)
        #return(ComplexHeatmap::draw(plot, merge_legends = TRUE, legend_grouping = "original"))
    }
}


















# Function to recover the list of genes found in pathways enriched
# @enrichResultObject: enrichResult object
# @column: column to parse
ParseEnrichResults <- function(enrichResultObject,
                               column = "geneID"#,
                               #output = "list"
                               )
{
    require("tidyr")
    if (!is.data.frame(enrichResultObject))
    {
        enrichResultObject <- data.frame(enrichResultObject)
    }
    #if (!output %in% c("list", "dataframe"))
    #{
    #    stop("The output argument must be either 'list' (the default) or 'dataframe'")
    #}
    #if (output == "list")
    #{
    #    l <- c()
    #    for (i in 1:nrow(enrichResultObject)){
    #        l <- c(l, str_split(enrichResultObject[[column]][i], "/")[[1]])
    #    }
    #    return(l[!duplicated(l)])
    #}
    #if (output == "dataframe")
    #{
   
    dataframe <- separate_rows(enrichResultObject,
                               core_enrichment,
                               sep = "/")
    dataframe <- data.frame(dataframe[c("Description", "core_enrichment")])
    return(dataframe)    
    #}
}


# Function from the mclust package
# ALlows to quantify how good is a clustering
adjustedRandIndex <- function (x, y)
    {
        x <- as.vector(x)
        y <- as.vector(y)
        if (length(x) != length(y))
            stop("arguments must be vectors of the same length")
        tab <- table(x, y)
        if (all(dim(tab) == c(1, 1)))
            return(1)
        a <- sum(choose(tab, 2))
        b <- sum(choose(rowSums(tab), 2)) - a
        c <- sum(choose(colSums(tab), 2)) - a
        d <- choose(sum(tab), 2) - a - b - c
        ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b +
                                                             a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
        return(ARI)
    }




################
# PCA function #
################

pca_func <- function(matrix,
                     sampleTable                 = NULL,
                     title                       = "PCA",
                     highlightedVar              = NULL,
                     pointSize                   = 8,
                     fontSize                    = 10,
                     highlightedVarName          = highlightedVar,
                     colors                      = NULL,
                     sampleLabel                 = NULL,
                     sampleLabelColumnCondition  = NULL,
                     sampleLabelValueCondition   = NULL,
                     sampleLabelSize             = 10,
                     stat_ellipse = FALSE
                     )
{
    # Load packages
    require(factoextra)
    require(ggrepel)
    # PCA analysis
    PC        <- prcomp(x = data.frame(t(matrix)),
                        center = TRUE,
                        scale. = TRUE
                        )
    eigVal    <- get_eigenvalue(PC)
    eigValPer <- round(eigVal$variance.percent, digits = 2)

    # Combine PCA result with sample plan (if sample plan is provided)
    if (!is.null(sampleTable))
    {
        pci <- data.frame(PC$x,sampleTable)
    } else {
        pci <- data.frame(PC$x)
    }
        
    # Plot PCA
    g <- ggplot(data = pci,
           mapping = aes(x = PC1, y = PC2)
           ) +
        xlab(paste0("PC1 (",eigValPer[1]," %)")) +
        ylab(paste0("PC2 (",eigValPer[2]," %)")) +
        ggtitle(title) +
        theme_bw()+
        theme(text              = element_text(size = fontSize),
              legend.key.height = unit(1.5,"cm")
              )
    if (!is.null(sampleLabel))
    {
        if (is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
        {
            g <- g +
                geom_text_repel(mapping = aes(label = eval(parse(text=sampleLabel))),
                                size = sampleLabelSize
                                )
        }
        else if (!is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
        {
            stop("missing sampleLabelValueCondition, it is required since sampleLabelColumnCondition is provided", call.=TRUE)
        }
        else if (!is.null(sampleLabelColumnCondition) && !is.null(sampleLabelValueCondition))
        {
            if (!is.element(el = sampleLabelValueCondition, set = as.vector(pci[,sampleLabelColumnCondition])))
            {
                warning(paste(sampleLabelColumnCondition,"column does not contain",sampleLabelValueCondition,"value. Therefore, no sample is labeled", sep = " "), call.=TRUE)
            }
            g <- g +
                geom_text_repel(data = dplyr::filter(.data = pci,
                                                     eval(parse(text=sampleLabelColumnCondition)) == sampleLabelValueCondition
                                                     ),
                                mapping = aes(label = eval(parse(text=sampleLabel))),
                                size = sampleLabelSize
                                )
        }
    }
    
    # Highlight variable
    if (!is.null(highlightedVarName))
    {
#        if (!is.element(highlightedVar,sampleTable))
#        {
#            stop("highlightedVar should be a valid column in sampleTable", call.=TRUE)
#        }
        g <- g + geom_point(aes(fill = eval(parse(text=highlightedVar))),
                            size = pointSize,
                            pch = 21
                            ) +
                                scale_fill_discrete(name = highlightedVar)
        if (!is.null(colors))
            {
                if (!is.null(names(colors)))
                {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors[levels(factor(umap_layout_plus[,highlightedVar]))]
                                               )
                } else {
                    g <- g + scale_fill_manual(name   = highlightedVarName,
                                               values = colors
                                               )
                }
            }
    } else {
        g <- g + geom_point(size = pointSize,
                            fill = "blue",
                            pch = 21
                            )
    }

    # Add ellipse
    if (stat_ellipse)
    {
        g <- g + stat_ellipse()
    }
    
    return(g)
}



#################
# tSNE function #
#################

tsne_func <- function(matrix,
                      sampleTable                = NULL,
                      highlightedVar             = NULL,
                      highlightedVarName         = highlightedVar,
                      colors                     = NULL,
                      seed                       = NULL,
                      perplexity                 = 30,
                      sampleLabelSize            = 4,
                      sampleLabel                = NULL,
                      sampleLabelColumnCondition = NULL,
                      sampleLabelValueCondition  = NULL,
                      fontsize                   = 30,
                      stat_ellipse               = FALSE
                      )
{
    # Load packages
    require(Rtsne)
    require(ggplot2)
    require(ggrepel)
    # Set seed
    set.seed(seed = seed)
    
    # Comput t-SNE
    tsne_model <- Rtsne(X                = t(matrix),
                        dims             = 2,
                        perplexity       = perplexity,
                        theta            = 0.0,
                        pca              = TRUE,
                        check_duplicates = FALSE,
                        is_distance      = FALSE
                        )

    # Get t-SNE matrix
    tsne_matrix <- as.data.frame(tsne_model$Y)

    # Bind sample plan
    if (!is.null(sampleTable))
        {
            tsne_matrix_plus <- cbind(sampleTable, tsne_matrix)
        } else {
            tsne_matrix_plus <- tsne_matrix
        }
    
    # Plot t-SNE
    g <- ggplot(data = tsne_matrix_plus,
                mapping = aes(x = V1, y = V2)
                ) +
                    geom_point(size = 6,
                               pch = 21,
                               fill = "blue"
                               ) +
                                   xlab("tSNE1") +
                                       ylab("tSNE2") +
                                           ggtitle("tSNE") +
                                               theme_bw()+
                                                   theme(text = element_text(size = fontsize),
                                                         legend.title = element_text(face = "bold"),
                                                         legend.key.height = unit(1.5,"cm")
                                                         )
    if (!is.null(sampleLabel))
        {
            if (is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
                {
                    g <- g +
                        geom_text_repel(mapping = aes(label = eval(parse(text=sampleLabel))),
                                        size = sampleLabelSize
                                        )
                }
            else if (!is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
                {
                    stop("missing sampleLabelValueCondition, it is required since sampleLabelColumnCondition is provided", call.=TRUE)
                }
            else if (!is.null(sampleLabelColumnCondition) && !is.null(sampleLabelValueCondition))
                {
                    if (!is.element(el = sampleLabelValueCondition, set = as.vector(tsne_matrix_plus[,sampleLabelColumnCondition])))
                        {
                            warning(paste(sampleLabelColumnCondition,"column does not contain",sampleLabelValueCondition,"value. Therefore, no sample is labeled", sep = " "), call.=TRUE)
                        }
                    g <- g +
                        geom_text_repel(data = dplyr::filter(.data = tsne_matrix_plus,
                                            eval(parse(text=sampleLabelColumnCondition)) == sampleLabelValueCondition
                                                             ),
                                        mapping = aes(label = eval(parse(text=sampleLabel))),
                                        size    = fontsize/5
                                        )
                }
        }

    # Highlight choosen variable
    if (!is.null(highlightedVar))
        {
            g <- g  + geom_point(aes(fill = eval(parse(text=highlightedVar))),
                                 size = 6,
                                 pch = 21
                                 ) + scale_fill_discrete(name = highlightedVar)
            if (!is.null(colors))
                {
                    if (!is.null(names(colors)))
                        {
                            g <- g + scale_fill_manual(name   = highlightedVarName,
                                                       values = colors[levels(factor(umap_layout_plus[,highlightedVar]))]
                                                       )
                        } else {
                            g <- g + scale_fill_manual(name   = highlightedVarName,
                                                       values = colors
                                                       )
                        }
                }
        } else {
            g <- g  +
                geom_point(size = 6,
                           fill = "blue",
                           pch = 21
                           )
        }
    
    # Add ellipse
    if (stat_ellipse)
        {
            g <- g + stat_ellipse()
        }    
    print(g)
    
}






plotByGene <- function(data,
                       Gene,
                       Annotations   = annotations,
                       group         = "Sample_Group",
                       main          = NULL,
                       whichMeasure  = "M values",
                       facet         = NULL)
    {
        require("ggplot2")
        require("ggrepel")

#        if (colnames(data) != rownames(annotations))
#            {
#                stop("please provide a data matrix with the same sample names as in the annotation (as rownames)")
#            }
        if (is.null(main))
        {
            main = Gene
        }

        if (is.null(facet)){            
            a <- data.frame(row.names    = rownames(Annotations),
                            Sample_Group = Annotations[[group]])
            a <- merge(a, t(data[Gene, ]), by = "row.names")
            names(a) <- c("Sample", "Sample_Group", "Beta")

        } else {
            print("I am here")
            a <- data.frame(row.names    = rownames(Annotations),
                        Sample_Group = Annotations[[group]],
                        facet        = Annotations[[facet]])
            a <- merge(a, t(data[Gene, ]), by = "row.names")
            names(a) <- c("Sample", "Sample_Group", facet, "Beta")

        }
        print(head(a))
        plot <- ggplot(a, aes(x=Sample_Group, y=Beta, fill=factor(Sample_Group))) +
            geom_boxplot(alpha=0.5) +
            #geom_jitter(alpha=0.2)
            #geom_point(position=position_jitterdodge(jitter.width=0.8), alpha=0.3) +
            stat_summary(fun=mean, colour="black", geom="point",
                         shape=18, size=3,show.legend = FALSE) +
            theme_bw() +
            ggtitle(main) +
            ylab(whichMeasure) + xlab("Sample group") + 
            labs(fill="Sample_Group") +
            theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle=45, hjust=1)) +
            geom_text_repel(aes(x=Sample_Group, y=Beta, label = Sample),
                            size = 2.5,
                            segment.color = "grey50",
                            position=position_jitterdodge(jitter.width=0.8))
        if (! is.null(facet))
        {
            plot <- plot + facet_wrap(facet)
        }
            #geom_text_repel(
            #    #data = subset(data, data$Sample_Group == "ATM" & data$Beta > mean_ATM),
            #    aes(label = Sample),
            #    size = 2.5,
            #    segment.color = "grey50",
            #    box.padding = unit(0.35, "lines"),
            #    point.padding = unit(0.6, "lines"),
            #    position=position_jitterdodge()       
        return(plot)
    }


write.gmt <- function(file,
                      type,
                      level,
                      sep,
                      colToUse = NULL)
{
    if (is.null(colToUse))
    {
        if (type == "GSEA")
        {
            colToUse <- "core_enrichment"
        } else {
            colToUse <- "geneID"
        }
    }
    obj <- read.csv(file, sep=sep)
    print(head(obj))
#    return(obj)
    if (level == "pathway")
    {
        for (analysis in unique(obj$comparison))
        {
            dataSubset <- subset(obj, obj$comparison == analysis)
            if (type == "GSEA"){
                NESSubset <- dataSubset[,c("Description", "ID", "NES")]
                dataSubset <- dataSubset[,c("Description", 'ID', colToUse, "NES")]
            } else {
                dataSubset <- dataSubset[,c("Description", 'ID', colToUse)]
            }
            dataSubset[[colToUse]] <- stringr::str_replace_all(dataSubset[[colToUse]], "/", "\t")
            print(head(dataSubset))

            if (type == "GSEA"){
                write.table(NESSubset,
                        file.path(dirname(file),
                                  paste0(stringr::str_split(basename(file), ".csv")[[1]][1],
                                         "_Pathways_",
                                         analysis,
                                         "_NES.gmt")
                                  ),
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE)
            }
         
            write.table(dataSubset,
                        file.path(dirname(file),
                                  paste0(stringr::str_split(basename(file), ".csv")[[1]][1],
                                         "_Pathways_",
                                         analysis,
                                         ".gmt")
                                  ),
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE)
        }
    }
    if (level == "gene")
    {
        for (analysis in unique(obj$comparison))
        {
            dataSubset <- subset(obj, obj$comparison == analysis)
            if (type == "GSEA"){
                dataSubset <- dataSubset[,c("Description", 'ID', colToUse, "NES")]
            } else {
                dataSubset <- dataSubset[,c("Description", 'ID', colToUse)]
            }
            dataSubset <- tidyr::separate_rows(dataSubset,
                                               all_of(colToUse), sep = "/")
            dataSubset <- data.frame(dataSubset)
            print(head(dataSubset))
            print(paste0(stringr::str_split(basename(file), ".csv")[[1]][1],
                                   "_Genes_",
                                   analysis,
                                   ".gmt")
                            )
            print(file.path(dirname(file),
                            paste0(stringr::str_split(basename(file), ".csv")[[1]][1],
                                   "_Genes_",
                                   analysis,
                                   ".gmt")
                            )
                  )
            write.table(dataSubset,
                        file.path(dirname(file),
                                  paste0(stringr::str_split(basename(file), ".csv")[[1]][1],
                                         "_Genes_",
                                         analysis,
                                         ".gmt")
                                  ),
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE)
        }
    }
}
