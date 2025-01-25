#!/usr/bin/env Rscript
#Torque Configuration
#PBS -l walltime=20:10:30
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=2
#PBS -q batch
#PBS -N Interpretation
#PBS -j oe
#PBS -o /data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation/


#---- Script du 29.10.2021 ----#
#---- Script mis à jour afin de réaliser les analyses d'enrichissement sur les nouveaux fichiers sans le duplicat (après normalisation Functional) ----#


# Trouver le dossier contenant le résultat des différentes comparaisons
DirDAGW <- file.path(Origin, "differentialAnalysis/GenomeWide/")
# Fichier d'annotation des gènes et promoteurs
DirMatrices = file.path(Origin, "matrices/")


th.p.val <- 0.05
th.logFC <- 1

# Fichier de destination de tous les résutats
#OutFiles1 = "~/ATM_Analysis/svn/analyse/results/FinalMethylation/Interpretation"
OutFiles1 = file.path(Origin, "Interpretation/")

# ACtivate renv project
renv::activate(project = "~/ATM_Analysis/svn/analyse/script/")

if(!dir.exists(OutFiles1))
{
    dir.create(OutFiles1, recursive = TRUE)
}


# Chargement des packages
library("stringr", quietly = TRUE)
library("pheatmap", quietly = TRUE)
library("ggplot2", quietly = TRUE)
library("viridis", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("tidyr", quietly = TRUE)
library("clusterProfiler", quietly = TRUE)
library("org.Hs.eg.db", quietly = TRUE)
library("enrichplot", quietly = TRUE)

setwd(DirDAGW)
files <- list.files(path = ".", pattern="_DMPs.csv", recursive = FALSE)


#######################################################
#################   Data retrieval  ###################
#######################################################

# Read the annotation file
annotations <- read.csv(file.path(DirMatrices, "Annotations.csv"), row.names = 1)

source(file = file.path(function.path, "functions.R"))


# create data frame with 0 rows and 5 columns
data <- data.frame(matrix(ncol = 7, nrow = 0))
# provide column names
colnames(data) <- c('Gene', 'logFC', 'P.Value', 'P.Value', 'adj.P.Val', 'comparison', 'Feature')
# Fill the dataframe
for (file in files){
    if (strsplit(file, "_")[[1]][2] == "Genes" | strsplit(file, "_")[[1]][2] == "Promoters"){
        print(file)
        a <- read.csv(file)
        #print(head(a))
        toAdd <- data.frame(Gene = a$X, logFC = a$logFC, P.Value = a$P.Value, adj.P.Val = a$adj.P.Val, comparison = strsplit(file, "_")[[1]][1], Feature = strsplit(file, "_")[[1]][2])
        data <- rbind(data, toAdd)
    }
}


#library("pathfindR")

#test <- input_processing(data2 %>%
#                         dplyr::filter(comparison == "A4" & Feature == "Promoters") %>%
#                         dplyr::select(Gene, logFC, adj.P.Val),
#                         pin_name_path="STRING")


#test

#output_df <- run_pathfindR(data2 %>%
#                           dplyr::filter(comparison %in% c("A1","A2","A3", "A4") & Feature == "Promoters") %>%
#                           dplyr::select(Gene, logFC, adj.P.Val),
#                           output_dir = "~/test_pathfindR",
#                           iterations = 20)



#test <- annotation <- biocarta_list <- fetch_gene_set(gene_sets = "BioCarta",
#                                min_gset_size = 10,
#                                max_gset_size = 300)

#biocarta_gsets <- biocarta_list[[1]]
#biocarta_descriptions <- biocarta_list[[2]]


#n_iter <- 10 ## number of iterations
#combined_res <- NULL ## to store the result of each iteration

#for (i in 1:n_iter)
#{
###### Active Subnetwork Search
#    snws_file <- paste0("active_snws_", i) # Name of output file
#    active_snws <- active_snw_search(input_for_search = RA_processed,
#                                     pin_name_path = "Biogrid",
#                                     snws_file = snws_file,
#                                     score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
#                                     sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
#                                     search_method = "GR")
###### Enrichment Analyses
#    current_res <- enrichment_analyses(snws = active_snws,
#                                       sig_genes_vec = RA_processed$GENE,
#                                       pin_name_path = "Biogrid",
#                                       genes_by_term = biocarta_gsets,
#                                       term_descriptions = biocarta_descriptions,
#                                       adj_method = "bonferroni",
#                                       enrichment_threshold = 0.05,
#                                       list_active_snw_genes = TRUE) # listing the non-input active snw genes in output
###### Combine results via `rbind`
#    combined_res <- rbind(combined_res, current_res)
#}

## Filter data
# For genes
data.genes.A0 <- data %>%
    filter(comparison == "A0" & Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.genes.A0A2 <- data %>%
    filter(comparison == "A0A2" & Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.genes.A1 <- data %>%
    filter(comparison == "A1" & Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.genes.A2 <- data %>%
    filter(comparison == "A2" & Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.genes.A3 <- data %>%
    filter(comparison == "A3" & Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.genes.A4 <- data %>%
    filter(comparison == "A4" & Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.genes.AVUS <- data %>%
    filter(comparison == "AVUS" & Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.genes.AVUSA2 <- data %>%
    filter(comparison == "AVUSA2" & Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.genes.A2Inactive <- data %>%
    filter(comparison == "A2Inactive" & Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.genes.All <- data %>%
    filter(Feature == "Genes") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")

# For promoters
data.pro.A0 <- data %>%
    filter(comparison == "A0" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.pro.A0A2 <- data %>%
    filter(comparison == "A0A2" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.pro.A1 <- data %>%
    filter(comparison == "A1" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.pro.A2 <- data %>%
    filter(comparison == "A2" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.pro.A3 <- data %>%
    filter(comparison == "A3" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.pro.A4 <- data %>%
    filter(comparison == "A4" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.pro.All <- data %>%
    filter(Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")

data.pro.AOld <- data %>%
    filter(comparison == "AOld" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")

data.pro.AOldER <- data %>%
    filter(comparison == "AOldER" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")

data.pro.AVUS <- data %>%
    filter(comparison == "AVUS" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.pro.AVUSA2 <- data %>%
    filter(comparison == "AVUSA2" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")
data.pro.A2Inactive <- data %>%
    filter(comparison == "A2Inactive" & Feature == "Promoters") %>%
    dplyr::select("Gene", "logFC", "P.Value", "adj.P.Val")


# Filter data according to logFC and adjusted p-value
data2 <- data
data2$logFC <-  as.numeric(data2$logFC)
data2 <- subset(data2, data2$adj.P.Val <= th.p.val & abs(data2$logFC) >= th.logFC)

# Recover genes
genes.A0 <- subset(data2, data2$comparison == "A0" & data2$Feature == "Genes")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
genes.A0A2 <- subset(data2, data2$comparison == "A0A2" & data2$Feature == "Genes")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
genes.A1 <- subset(data2, data2$comparison == "A1" & data2$Feature == "Genes")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
genes.A2 <- subset(data2, data2$comparison == "A2" & data2$Feature == "Genes")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
genes.A3 <- subset(data2, data2$comparison == "A3" & data2$Feature == "Genes")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
genes.A4 <- subset(data2, data2$comparison == "A4" & data2$Feature == "Genes")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
genes.AVUS <- subset(data2, data2$comparison == "AVUS" & data2$Feature == "Genes")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
genes.AVUSA2 <- subset(data2, data2$comparison == "AVUSA2" & data2$Feature == "Genes")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
genes.A2Inactive <- subset(data2, data2$comparison == "A2Inactive" & data2$Feature == "Genes")[c("Gene", "logFC", "P.Value", "adj.P.Val")]


# Recover promoters
pro.A0 <- subset(data2, data2$comparison == "A0" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
pro.A0A2 <- subset(data2, data2$comparison == "A0A2" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
pro.A1 <- subset(data2, data2$comparison == "A1" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
pro.A2 <- subset(data2, data2$comparison == "A2" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
pro.A3 <- subset(data2, data2$comparison == "A3" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
pro.A4 <- subset(data2, data2$comparison == "A4" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]

pro.AOld <- subset(data2, data2$comparison == "AOld" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
pro.AVUS <- subset(data2, data2$comparison == "AVUS" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
pro.AVUSA2 <- subset(data2, data2$comparison == "AVUSA2" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]
pro.A2Inactive <- subset(data2, data2$comparison == "A2Inactive" & data2$Feature == "Promoters")[c("Gene", "logFC", "P.Value", "adj.P.Val")]

# Recover genes and promoters
genesAndPro.A0 <- rbind(genes.A0, pro.A0)
genesAndPro.A0A2 <- rbind(genes.A0A2, pro.A0A2)
genesAndPro.A1 <- rbind(genes.A1, pro.A1)
genesAndPro.A2 <- rbind(genes.A2, pro.A2)
genesAndPro.A3 <- rbind(genes.A3, pro.A3)
genesAndPro.A4 <- rbind(genes.A4, pro.A4)
genesAndPro.AVUS <- rbind(genes.AVUS, pro.AVUS)
genesAndPro.AVUSA2 <- rbind(genes.AVUSA2, pro.AVUSA2)
genesAndPro.A2Inactive <- rbind(genes.A2Inactive, pro.A2Inactive)

res1 <- data.frame(row.names = unique(data2$Gene),
           "A1" = unique(data2$Gene) %in% c(pro.A1$Gene, genes.A1$Gene),
           "A2" = unique(data2$Gene) %in% c(pro.A2$Gene, genes.A2$Gene),
           "A3" = unique(data2$Gene) %in% c(pro.A3$Gene, genes.A3$Gene),
           "A4" = unique(data2$Gene) %in% c(pro.A4$Gene, genes.A4$Gene),
           "A0" = unique(data2$Gene) %in% c(pro.A0$Gene, genes.A0$Gene),
           "A0A2" = unique(data2$Gene) %in% c(pro.A0A2$Gene, genes.A0A2$Gene),
           "AVUS" = unique(data2$Gene) %in% c(pro.AVUS$Gene, genes.AVUS$Gene),
           "AVUSA2" = unique(data2$Gene) %in% c(pro.AVUSA2$Gene, genes.AVUSA2$Gene),
           "A2Inactive" = unique(data2$Gene) %in% c(pro.A2Inactive$Gene, genes.A2Inactive$Gene))

res1[res1 == "TRUE"] <- 1
res1[res1 == 0] <- 1

res2 <- data.frame(
    row.names = unique(data2$Gene),
    Feature = ifelse(unique(data2$Gene) %in% c(pro.A1$Gene, pro.A2$Gene, pro.A3$Gene, pro.A4$Gene, pro.A0$Gene, pro.A2Inactive$Gene),
        ifelse(unique(data2$Gene) %in% c(genes.A1$Gene, genes.A2$Gene, genes.A3$Gene, genes.A4$Gene, genes.A0$Gene, genes.A2Inactive$Gene),
               "ProAndGene", "Promoters"), "Genes")
)


write.csv(res1, file.path(OutFiles1, "AnalysisInfo_cytoscape.csv"))
write.csv(res2, file.path(OutFiles1, "FeatureInfo_cytoscape.csv"))



########################################################
##############  Functions for GSEA plots  ##############
########################################################



# Plot the GSEA results
gsInfo <- function(object, geneSetID) {
        geneList <- object@geneList

            if (is.numeric(geneSetID))
                        geneSetID <- object@result[geneSetID, "ID"]

            geneSet <- object@geneSets[[geneSetID]]
            exponent <- object@params[["exponent"]]
            df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
            df$ymin <- 0
            df$ymax <- 0
            pos <- df$position == 1
            h <- diff(range(df$runningScore))/20
            df$ymin[pos] <- -h
            df$ymax[pos] <- h
            df$geneList <- geneList

            df$Description <- object@result[geneSetID, "Description"]
            return(df)
        }

gseaScores <- getFromNamespace("gseaScores", "DOSE")


tableGrob2 <- function(d, p = NULL) {
        # has_package("gridExtra")
        d <- d[order(rownames(d)),]
        tp <- gridExtra::tableGrob(d)
        if (is.null(p)) {
                    return(tp)
                        }

        # Fix bug: The 'group' order of lines and dots/path is different
        p_data <- ggplot_build(p)$data[[1]]
        # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
        p_data <- p_data[order(p_data[["group"]]), ]
        pcol <- unique(p_data[["colour"]])
        ## This is fine too
        ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]
        j <- which(tp$layout$name == "rowhead-fg")

        for (i in seq_along(pcol)) {
                    tp$grobs[j][[i+1]][["gp"]] <- grid::gpar(col = pcol[i])
                        }
        return(tp)
    }

gseaplot2 <- function (x, geneSetID, title = "", color = "green", base_size = 11, 
                       rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
                       ES_geom = "line") 
{
    require(RColorBrewer)
    ES_geom <- match.arg(ES_geom, c("line", "dot"))
    geneList <- position <- NULL
    if (length(geneSetID) == 1) {
        gsdata <- gsInfo(x, geneSetID)
    }
    else {
        gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
    }
    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
        theme(panel.grid.major = element_line(colour = "grey92"), 
            panel.grid.minor = element_line(colour = "grey92"), 
            panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
        scale_x_continuous(expand = c(0, 0))
    if (ES_geom == "line") {
        es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
            size = 1)
    }
    else {
        es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
            size = 1, data = subset(gsdata, position == 1))
    }
    p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
        legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
    p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
            unit = "cm"))
    i <- 0
    for (term in unique(gsdata$Description)) {
        idx <- which(gsdata$ymin != 0 & gsdata$Description == 
            term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
    }
    p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
        ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
        theme_classic(base_size) + theme(legend.position = "none", 
        plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
        axis.text = element_blank(), axis.line.x = element_blank()) + 
        scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
        0))
    if (length(geneSetID) == 1) {
        v <- seq(1, sum(gsdata$position), length.out = 9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0) 
            inv <- inv + 1
        col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
        ymin <- min(p2$data$ymin)
        yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
            xmax = xmax, col = col[unique(inv)])
        p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
            ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
            alpha = 0.9, inherit.aes = FALSE)
    }
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
        y = ~y, yend = 0), color = "grey")
    p.pos <- p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") + 
        theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
            l = 0.2, unit = "cm"))
    if (!is.null(title) && !is.na(title) && title != "") 
        p.res <- p.res + ggtitle(title)
    if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values = color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        }
        else {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }
    if (pvalue_table) {
        pd <- x[geneSetID, c("Description", "enrichmentScore", "NES", "pvalue", "p.adjust")]
        rownames(pd) <- pd$Description
        names(pd)[names(pd) == "enrichmentScore"] <- "ES"
        pd <- pd[, -1]
        #pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
        #rownames(pd) <- pd$Description
        #pd <- pd[, -1]
        pd <- round(pd, 4)
        tp <- tableGrob2(pd, p.res)
        p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
            xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
                0.95), ymin = quantile(p.res$data$runningScore, 
                0.75), ymax = quantile(p.res$data$runningScore, 
                0.9))
    }
    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
        axis.ticks.x = element_line(), axis.text.x = element_text())
    if (length(subplots) == 1) 
        return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
            r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
    if (length(rel_heights) > length(subplots)) 
        rel_heights <- rel_heights[subplots]
    aplot::plot_list(gglist = plotlist, ncol = 1, heights = rel_heights)
}

# plotGSEA: function to plot GSEA results but also a heatmap and umap of ecah enriched pathway
# @ListOfEnrichmentObject: a list of names of enrichment objects to use (indicated as string)
plotGSEA <- function(ListOfEnrichmentObject){
    for (i in 1:length(ListOfEnrichmentObject))
    {
        # Get informations for the names of the outputs
        Feature <- print(strsplit(ListOfEnrichmentObject[i], "\\.")[[1]][2])
        Analysis <- print(strsplit(ListOfEnrichmentObject[i], "\\.")[[1]][1])
        # Evaluate the object (from name to GseaResult)
        EnrichmentObject <- print(eval(parse(text = ListOfEnrichmentObject[i])))
        name <- paste0(Analysis, ".", Feature)
        print(name)
        
        for (pathN in 1:nrow(EnrichmentObject)){
            pathway <- stringr::str_replace(EnrichmentObject$Description[pathN],
                                            "/",
                                            " ")
            # GSEA plot 
            gseaN <- gseaplot2(EnrichmentObject,
                               geneSetID     = pathN,
                               title         = EnrichmentObject$Description[pathN],
                               pvalue_table  = TRUE)
            ggsave(file.path(OutFiles,
                             paste0(name, "_GseaPlot_",
                                    pathway,
                                    '.png')),
                   gseaN,
                   height  = 7,
                   width   = 14)
            # Get the list of genes for the enriched terms
            List <- ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description == EnrichmentObject$Description[pathN],],
                                       column = "core_enrichment")
            List <- List$core_enrichment
            
            # Heatmap of corresponfing genes
            if (Feature == "Pro")
            {
                dataFrame.annotated = M.Promoters.annotated
                dataFrame.annotated2 = Beta.Promoters.annotated
            }
            if (Feature == "Genes")
            {
                dataFrame.annotated = M.Genes.annotated
                dataFrame.annotated2 = Beta.Genes.annotated
            }            
                
            try(
                heatmapN <- heatmap_func(data          = dataFrame.annotated[List,],
                                         annotations   = annotations,
                                         fontsize_row  = 8,
                                         fontsize_col  = 10,
                                         main          = EnrichmentObject$Description[pathN])
            )
            try(
                ggsave(file.path(OutFiles,
                                 paste0(name, "_M_Heatmap_",
                                        pathway,
                                        '.png')),
                       heatmapN,
                       height  = 7,
                       width   = 14)
            )
            
            try(
                heatmapN2 <- heatmap_func(data          = dataFrame.annotated2[List,],
                                     annotations   = annotations,
                                     fontsize_row  = 8,
                                     fontsize_col  = 10,
                                     main          = EnrichmentObject$Description[pathN])
            )
            try(
                ggsave(file.path(OutFiles,
                                 paste0(name, "_Beta_Heatmap_",
                                        pathway, '.png')),
                       heatmapN2,
                       height  = 7,
                       width   = 14)
            )
            # UMAP of corresponding genes
            k <- 2
            try(
                umap.plot <- umap_func(matrix  = dataFrame.annotated[List,],
                                       sampleTable         = annotations,
                                       highlightedVar      = "Sample_Group",
                                       highlightedVarName  = "Sample_Group",
                                       colors              = c("orange", "Red"),
                                       sampleLabel         = "Sample_Name",
                                       metric              = "euclidean",
                                       title               = EnrichmentObject$Description[pathN],
                                       fontsize            = 14,
                                       assignGroup         = TRUE,
                                       k = k
                                       )
            )
            try(
                ggsave(filename = file.path(OutFiles,
                                            paste0(name, "_UMAP_",
                                                   pathway,
                                                   '.png')),
                       plot = umap.plot,
                       width = 12, height = 7)
            )
        }
    }
}













##########################################################
#### Enrichemnt with c6.Oncogenic signature gene sets ####
##########################################################

# ~/ATM_Analysis/data/c6.all.v7.5.1.symbols.gmt

