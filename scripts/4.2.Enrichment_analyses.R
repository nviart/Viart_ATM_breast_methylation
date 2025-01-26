#!/usr/bin/env Rscript
#Torque Configuration
#PBS -l walltime=20:10:30
#PBS -l mem=45gb
#PBS -l nodes=2:ppn=4
#PBS -q batch
#PBS -N Interpretation_Kegg_universe
#PBS -j oe


# Source the file that will find data to use
source("~/ATM_Analysis/svn/analyse/script/methylation/MethPipeline/config.R")
renv::activate(project = renv.path)
source("~/ATM_Analysis/svn/analyse/script/methylation/MethPipeline/4.Interpretation.R")

#load(file.path(Origin, "/matrices/Beta_Funnorm_Normalisation.RData"))
#load(file.path(Origin, "/matrices/M_Funnorm_Normalisation.RData"))
load(file.path(Origin, "/matrices/Mapping.RData"))

# Get locatio of genes differentially methylated
cytoBands <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = keys(org.Hs.eg.db),
                                   columns = c("SYMBOL", "MAP"))

cytoSelected <- subset(cytoBands, cytoBands$SYMBOL %in% c(pro.A1$Gene, pro.A2$Gene))
cytoSelected  <- cytoSelected[order(cytoSelected$MAP), ]

cytoSelected$A1 <- ifelse(cytoSelected$SYMBOL %in% pro.A1$Gene, "True", "False")
cytoSelected$A2 <- ifelse(cytoSelected$SYMBOL %in% pro.A2$Gene, "True", "False")

write.csv(cytoSelected, file.path(Origin, "cytogeneticBands.csv"))

# Add annotations of LOH (from Anne-Laure paper)
ATM.LOH.Yes <- c("T0220", "T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
ATM.LOH.No <- c("T0120", "T0015", "T0009", "T0249", "T0247")

annotations$ATM_LOH  <- ifelse(annotations$Sample_Name %in% ATM.LOH.Yes, "LOH",
                               ifelse(annotations$Sample_Name %in% ATM.LOH.No, "No LOH", "NA")
                               )


# Load the universe genes and promoters
geneUniverse <- read.csv(file.path(Origin, "matrices/probes_Genes.csv"), row.names=1)
geneUniverse <- unique(geneUniverse$SYMBOL)

proUniverse <- read.csv(file.path(Origin, "matrices/probes_Promoters.csv"), row.names=1)
proUniverse <- unique(proUniverse$SYMBOL)



max <- max(c(length(pro.A0A2$Gene), length(pro.A2$Gene), length(pro.AVUSA2$Gene), length(pro.A0$Gene), length(pro.A1$Gene)))

allData <- data.frame(
    "ATM PV ER+" = c(pro.A2$Gene, rep("", max - length(pro.A2$Gene))),
    "ATM Inactive ER+" = c(pro.A2Inactive$Gene, rep("", max - length(pro.A2Inactive$Gene))),
    "All ER+" = c(pro.A0A2$Gene, rep("", max - length(pro.A0A2$Gene))),
    "ATM VUS ER+" = c(pro.AVUSA2$Gene, rep("", max - length(pro.AVUSA2$Gene))),
    A0 = c(pro.A0$Gene, rep("", max - length(pro.A0$Gene))),
    A1 = c(pro.A1$Gene, rep("", max - length(pro.A1$Gene)))
    )


write.csv(allData, file.path(Origin, "differentialAnalysis", "DM_promoters.csv"))





oncoKB <- read.csv("/data/kdi_prod/.kdi/project_workspace_0/833/acl/05.00/data/establishedCancerGenes/OncoKB_cancerGeneList_update2024.03.21.tsv", sep="\t")
oncoKB <- c(oncoKB[["Hugo.Symbol"]],
            data.frame(tidyr::separate_rows(oncoKB,
                                            Gene.Aliases,
                                            sep = ","))[["Gene.Aliases"]]
            )
oncoKB <- oncoKB[!duplicated(oncoKB)]
oncoKB <- stringr::str_replace(oncoKB, " ", "")


pro.A2$Gene[pro.A2$Gene %in% oncoKB]

pro.A2Inactive$Gene[pro.A2Inactive$Gene %in% oncoKB]

head(oncoKB)



dim(pro.A2Inactive[pro.A2Inactive$logFC > 0,])

dim(pro.A2Inactive[pro.A2Inactive$logFC < 0,])


dim(pro.A0A2[pro.A0A2$logFC > 0,])

dim(pro.A0A2[pro.A0A2$logFC < 0,])





#######################################################
#################   MSIG database  ####################
#######################################################

OutFiles = file.path(OutFiles1, "MSIGDB/Gse/")

if(!dir.exists(OutFiles))
{
    dir.create(OutFiles, recursive = TRUE)
}


getMsigdbGSEA <- function(geneList,
                          comparison = NULL,
                          map        = FALSE,
                          seed       = NULL,
                          NbrProbes  = NULL)
{
    require("msigdbr")

    returnENTREZID <- function(gse, corr){
        for (i in 1:dim(gse)[1])
        {
            toMap <- str_split(gse$core_enrichment[i], "\\s*/\\s*")[[1]]
            mapped <- corr[corr$ENTREZID %in% toMap, ]$SYMBOL
            #print(length(toMap) == length(mapped))
            gse$core_enrichment_SYMBOL[i] <- paste(mapped, collapse = "/")
        }
        return(gse)
    }

    # Set seed
    set.seed(1234)
    feature <- gsub(".*\\.(.*)\\..*", "\\1", deparse(substitute(geneList)))

    if (!is.null(NbrProbes) && feature == "pro")
    {
        # Retrieve the number of probes in promoters and genes
        NbrProbesPromoters <- read.csv(file.path(Origin, "matrices/M_Promoterss_NbrProbes.csv"),
                                       row.names = 1)
        print("Number of probes read")
        geneList <- merge(geneList, NbrProbesPromoters, by.x="Gene", by.y="row.names")
        geneList$logFC <- geneList$NbrProbe * geneList$logFC
    }
    
    # Query Msigdb
    msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
        dplyr::select(gs_name, entrez_gene) %>%
            dplyr::rename(ont = gs_name, gene = entrez_gene)

    print("Number of gene sets tested: ") 
    print(length(unique(msig_h$ont)))

    # Create the dataframe needed for GSEA analysis
    corr <- bitr(geneList$Gene,
                 fromType  = "SYMBOL",
                 toType    = "ENTREZID",
                 OrgDb     = "org.Hs.eg.db")

    d <- merge(corr, geneList, by.x = "SYMBOL", by.y = "Gene")[c("ENTREZID", "logFC")]
    geneList = d[,2]
    names(geneList) = as.character(d[,1])
    geneList = sort(geneList, decreasing = TRUE)

#    d <- geneList
#    geneList <- d$logFC
#    names(geneList) <- as.character(d$Gene)
#    geneList = sort(geneList, decreasing = TRUE)
#    print(head(geneList))

    # Perform GSEA analysis
    tryCatch(
        expr = {
            Gse_MSIG <- GSEA(gene          = geneList,
                             TERM2GENE      = msig_h,
                             minGSSize      = 2,
                             maxGSSize      = 500,
                             pvalueCutoff   = 0.05,
                             pAdjustMethod  = "BH",
                             verbose        = TRUE,
                             seed           = TRUE) %>%
                as_tibble
        
    
    #Gse_MSIG <- setReadable(Gse_MSIG,
    #                        OrgDb   = org.Hs.eg.db,
    #                        keyType = "ENTREZID")
        
            Gse_MSIG <- Gse_MSIG %>%
                arrange(desc(abs(NES)))
            Gse_MSIG <- returnENTREZID(Gse_MSIG, corr)
            return(Gse_MSIG)
        },
        error = function(e){
            return(NULL)
        }
    )
}



#msig_comp <- lolliplot_custom_Paper(data = rbind(as.data.frame(A1.Pro.msig),
#                                                 as.data.frame(A2.Pro.msig)),
#                                    comparison = TRUE,
#                                    main              = NULL,
#                                    orderAccordingTo  = "GeneRatio",
#                                    size              = "NES",
#                                    showSign          = T,
#                                    showCategory      = NULL,
#                                    sign              = NULL
#                                    )

A0A2.Pro.msig <- getMsigdbGSEA(data.pro.A0A2,
                               comparison = NULL,
                               map        = FALSE,
                               seed       = NULL)
A0A2.Pro.msig$comparison <- "A0A2"


A2.Pro.msig <- getMsigdbGSEA(data.pro.A2,
                             comparison = NULL,
                             map        = FALSE,
                             seed       = NULL)
A2.Pro.msig$comparison <- "A2"


A2Inactive.Pro.msig <- getMsigdbGSEA(data.pro.A2Inactive,
                             comparison = NULL,
                             map        = FALSE,
                             seed       = NULL)
A2Inactive.Pro.msig$comparison <- "A2Inactive"


write.csv(rbind(A2.Pro.msig, A2Inactive.Pro.msig, A0A2.Pro.msig),
          file.path(OutFiles, "msig_pathways.csv"))



A0A2.Pro.msig <- getMsigdbGSEA(data.pro.A0A2,
                               comparison = NULL,
                               map        = FALSE,
                               seed       = NULL,
                               NbrProbes  = TRUE)
A0A2.Pro.msig$comparison <- "A0A2"


A2.Pro.msig <- getMsigdbGSEA(data.pro.A2,
                             comparison = NULL,
                             map        = FALSE,
                             seed       = NULL,
                             NbrProbes  = TRUE)
A2.Pro.msig$comparison <- "A2"


A2Inactive.Pro.msig <- getMsigdbGSEA(data.pro.A2Inactive,
                             comparison = NULL,
                             map        = FALSE,
                             seed       = NULL,
                             NbrProbes  = TRUE)
A2Inactive.Pro.msig$comparison <- "A2Inactive"


write.csv(rbind(A2.Pro.msig, A2Inactive.Pro.msig, A0A2.Pro.msig),
          file.path(OutFiles, "Nbr_msig_pathways.csv"))

file <- file.path(OutFiles, "/Nbr_msig_pathways.csv")
write.gmt(file, type = "GSEA", level = "pathway", sep=",", colToUse="core_enrichment_SYMBOL")



#path <- c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_MITOTIC_SPINDLE",
#          "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
#          "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS",
#          "HALLMARK_NOTCH_SIGNALING", "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_UV_RESPONSE_UP",
#          "HALLMARK_GLYCOLYSIS", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
#          "HALLMARK_HEDGEHOG_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
#          "HALLMARK_MYC_TARGETS_V1", "HALLMARK_ADIPOGENESIS",
#          "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_HEME_METABOLISM",
#          "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
#          "HALLMARK_UV_RESPONSE_DN", "HALLMARK_APOPTOSIS", "HALLMARK_HYPOXIA",
#          "HALLMARK_ESTROGEN_RESPONSE_LATE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
#          "HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR", "HALLMARK_IL2_STAT5_SIGNALING",
#          "HALLMARK_PEROXISOME", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_COMPLEMENT",
#          "HALLMARK_PROTEIN_SECRETION", "HALLMARK_MYC_TARGETS_V2",
#          "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_APICAL_JUNCTION")

#cat <- c("DNA repair pathways", "Cancer related pathways", "Metabolism",
#         "Cellular growth/Division",
#         "Protein misfolding", "Immunology")

#label <- c(cat[4], cat[4],
#           cat[6], cat[3],
#           cat[4] , cat[4],
#           cat[2], cat[3], cat[1],
#           cat[3], cat[2],
#           cat[2], cat[5],
#           cat[2], cat[3],
#           cat[3], cat[3],
#           cat[6], cat[6],
#           cat[1], cat[2], cat[3],
#           cat[3], cat[2],
#           cat[2], cat[1], cat[6],
#           cat[3], cat[3], cat[6],
#           cat[3], cat[2],
#           cat[6], cat[2])
#names(label) <- path
#label <- data.frame(label)
#label <- cbind(label, "comp")
#names(label) <- c("label", "comparison")

#data = rbind(as.data.frame(A0.Pro.msig),
#             as.data.frame(A1.Pro.msig),
#             as.data.frame(A2.Pro.msig),
#             as.data.frame(AVUS.Pro.msig),
#             as.data.frame(AVUSA2.Pro.msig))



#ggplot(data.frame(label)) + geom_point(aes(x=1,
#                                           y=c(1:dim(data.frame(label))[1]),
#                                           colour = label, shape = label),
#                                       size=2)
           

comment <- function(){

msig_comp <- lolliplot_custom_Paper(data = rbind(as.data.frame(A0A2.Pro.msig),
                                                 as.data.frame(A2.Pro.msig),
                                                 as.data.frame(AVUSA2.Pro.msig)),
                                    #categories = label,
                                    comparison = TRUE,
                                    main              = NULL,
                                    orderAccordingTo  = "GeneRatio",
                                    size              = "NES",
                                    showSign          = T,
                                    showCategory      = NULL,
                                    sign              = NULL
                                    )

library(grid)
gt = ggplot_gtable(ggplot_build(msig_comp))

index <- gt$layout$l[grep('panel-1-1', gt$layout$name)]
# modify the width
gt$widths[index] = unit(0.01, "npc")#0.2*gt$widths[index]
# Suppress x ticks
index <- grep('axis-b-1-1', gt$layout$name)
gt$grobs[[index]]$children[2]$axis[1]$grob[[1]]$gp$col <- "NA"
gt$grobs[[index]]$children[2]$axis[1]$grob[[1]]$gp$fill <- "NA"
#<- nullGrob()
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




grid.draw(gt)


png(filename="~/test_lolliplot.png", width=1800, height=1000, units="px", res=100)
grid.draw(gt)
dev.off()

}




# Specific to A1 (ER+)
#match(A1.Pro.msig$Description, A2.Pro.msig$Description)
#A1.Pro.msig$Description[which(match(A1.Pro.msig$Description, A2.Pro.msig$Description) %in% NA)]

# Specific to A2 (ATM tumours)
#match(A2.Pro.msig$Description, A1.Pro.msig$Description)
#A2.Pro.msig$Description[which(match(A2.Pro.msig$Description, A1.Pro.msig$Description) %in% NA)]



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
        
        for (pathN in 1:nrow(EnrichmentObject)){
            # Get the list of genes for the enriched terms
            List <- ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description == EnrichmentObject$Description[pathN],],
                                       column = "core_enrichment_SYMBOL")
            List <- List$core_enrichment
            print(List)
            
            # Heatmap of corresponfing genes
            if (Feature == "Pro")
            {
                dataFrame.annotated = data_Promoters
                dataFrame.annotated2 = data_Genes
            }
            if (Feature == "Genes")
            {
                dataFrame.annotated = M.Genes.annotated
                dataFrame.annotated2 = Beta.Genes.annotated
            }            

            print("color selected")
            print(colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100))
            if ("clinvar" %in% colnames(annotations)){
                annotations[is.na(annotations$clinvar), ]$clinvar <- ""
            }
            heatmapN <- heatmap_func(data          = dataFrame.annotated[List,],
                                     annotations   = annotations,
                                     fontsize_row  = 8,
                                     fontsize_col  = 10,
                                     color         = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),

                                     main          = EnrichmentObject$Description[pathN],
                                     filename = file.path(OutFiles,
                                                          paste0('test.png')))
            print("first done...")

            heatmapN2 <- heatmap_func(data         = dataFrame.annotated2[List,],
                                     annotations   = annotations,
                                     fontsize_row  = 8,
                                     fontsize_col  = 10,
                                     color         = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),

                                     main          = EnrichmentObject$Description[pathN],
                                     filename = file.path(OutFiles,
                                                          paste0(name, "_Beta_Heatmap_",
                                                                 EnrichmentObject$Description[pathN], '.png')))
            
        }
    }
}


#plotGSEA("A0A2.Pro.msig")

#plotGSEA("A2.Pro.msig")











#######################################################
#################   KEGG database  ####################
#######################################################

# getKeggEnriched: retreive the enrichemnt of gene pathways
# @genes: a list of genes
# @comparison: the comparison to perform
getKeggEnriched <- function(genes,
                            comparison = NULL,
                            universe = NULL)
{
    # Create the dataframe needed
    corr <- bitr(genes,
                 fromType  = "SYMBOL",
                 toType    = "ENTREZID",
                 OrgDb     = "org.Hs.eg.db")

    if (!is.null(universe))
    {
        print("I am here")
        universe <- bitr(universe,
                         fromType  = "SYMBOL",
                         toType    = "ENTREZID",
                         OrgDb     = "org.Hs.eg.db")
        universe <- universe$ENTREZID
    }
    # Perform the analysis
    Enriched_Kegg <- enrichKEGG(gene           = corr$ENTREZID,
                                organism       = 'hsa',
                                keyType        = "ncbi-geneid",
                                #universe       = universe,
                                pvalueCutoff   = 0.05,
                                pAdjustMethod  = "BH",
                                qvalueCutoff   = 0.2,
                                minGSSize      = 2,
                                maxGSSize      = 500)
    print(Enriched_Kegg)
    # Return the dataframe
    if (nrow(data.frame(Enriched_Kegg)) != 0)
        {
            Enriched_Kegg <- setReadable(Enriched_Kegg,
                                         OrgDb = org.Hs.eg.db,
                                         keyType="ENTREZID")
            if (!is.null(comparison))
            {
                return(data.frame(Enriched_Kegg, comparison = comparison))
            }
            else
            {
                return(data.frame(Enriched_Kegg))
            }
        }
}



########################################
##### Over-representation analysis #####
########################################

# Results will be outputted to
OutFiles = file.path(OutFiles1, "Kegg/PathwayOver-representation/") 

if(!dir.exists(OutFiles))
{
    dir.create(OutFiles, recursive = TRUE)
}


#############   For Genes  ################

A2up.Genes.Kegg <- getKeggEnriched(genes.A2[genes.A2$logFC > 0,]$Gene, comparison = "A2")


# Retrieve Kegg for Genes
A2.Genes.Kegg <- getKeggEnriched(genes.A2$Gene, comparison = "A2")
A0A2.Genes.Kegg <- getKeggEnriched(genes.A0A2$Gene, comparison = "A0A2")
A2Inactive.Genes.Kegg <- getKeggEnriched(genes.A2Inactive$Gene, comparison = "A2Inactive")

# Plot the comparison of the 4 analysis
Kegg.all <- rbind(data.frame(A2.Genes.Kegg),
                  data.frame(A0A2.Genes.Kegg),
                  data.frame(A2Inactive.Genes.Kegg))
write.csv(Kegg.all, paste0(OutFiles, "Kegg_Genes_according_comparisons.csv"))

#if (nrow(Kegg.all) != 0){
#    plot.genes.Kegg <- dotplot_custom(Kegg.all, main = "Enriched Kegg pathways - Genes", comparison=TRUE)
#    ggsave(paste0(OutFiles, "Comparisons_Kegg_Genes.png"), plot.genes.Kegg, width= 16, height=8)
#}

# Plot for each level of enrichments with all genes (all comparisons confounded)
#if (is.null(all.Genes.Kegg) == FALSE){
#    write.csv(all.Genes.Kegg, paste0(OutFiles, "Kegg_ALLGenes.csv"))
#    plot.genes.Kegg.all <- dotplot_custom(all.Genes.Kegg, main = "Enriched Kegg pathways - all Genes")
#    ggsave(paste0(OutFiles, "Comparisons_Kegg_ALLGenes.png"), plot.genes.Kegg.all, width= 16, height=8)
#}

#rm(A2.Genes.Kegg, A3.Genes.Kegg, A4.Genes.Kegg, all.Genes.Kegg, Kegg.all)
#rm(Kegg.all, plot.genes.Kegg, plot.genes.Kegg.all)

#############   For Promoters  ################

## For Promoters
# Retrieve Kegg for Promoters
A2.Pro.Kegg <- getKeggEnriched(pro.A2$Gene, comparison = "A2")
A0A2.Pro.Kegg <- getKeggEnriched(pro.A0A2$Gene, comparison = "A0A2")
A2Inactive.Pro.Kegg <- getKeggEnriched(pro.A2Inactive$Gene, comparison = "A2Inactive")

# Plot the comparison of the 4 analysis
Kegg.Pro <- rbind(data.frame(A2.Pro.Kegg),
                  data.frame(A0A2.Pro.Kegg),
                  data.frame(A2Inactive.Pro.Kegg))
write.table(Kegg.Pro, paste0(OutFiles, "Kegg_Promoters_according_comparisons.csv"), sep="\t")

#if (nrow(Kegg.Pro) != 0){
#    plot.Pro.Kegg <- dotplot_custom(Kegg.Pro, main = "Enriched Kegg pathways - Promoters", comparison=TRUE)
#    ggsave(paste0(OutFiles, "Comparisons_Kegg_Pro.png"), plot.Pro.Kegg, width= 16, height=8)
#}

# Plot for each level of enrichments with all genes (all comparisons confounded)
#if (is.null(all.Pro.Kegg) == FALSE){
#    write.csv(all.Pro.Kegg, paste0(OutFiles, "Kegg_ALLPromoters.csv"))
#    plot.Pro.Kegg.all <- dotplot_custom(all.Pro.Kegg, main = "Enriched Kegg pathways - all Promoters")
#    ggsave(paste0(OutFiles, "Comparisons_Kegg_ALLPromoters.png"), plot.Pro.Kegg.all, width= 16, height=8)
#}

#rm(A1.Pro.Kegg, A2.Pro.Kegg, A3.Pro.Kegg, A4.Pro.Kegg, all.Pro.Kegg)
#rm(Kegg.Pro, plot.Pro.Kegg, plot.Pro.Kegg.all)












print("beginning of GSEA for KEGG")

##########################################
##### Gene set enrichments analysis  #####
##########################################

# Results will be outputted to
OutFiles = file.path(OutFiles1, "Kegg/Gse/") 

if(!dir.exists(OutFiles))
{
    dir.create(OutFiles, recursive = TRUE)
}


getKeggGSE <- function(geneList,
                       comparison    = NULL,
                       map           = FALSE,
                       seed          = NULL,
                       NbrProbes     = NULL,
                       pvalueCutoff  = 0.05)
{
    # Set seed
    #set.seed(seed)
    set.seed(1234)
    feature <- gsub(".*\\.(.*)\\..*", "\\1", deparse(substitute(geneList)))

    print(feature)

    if (!is.null(filter))
    {
        NbrProbes = FALSE
    }

    print(head(geneList))

    if (!is.null(NbrProbes) && feature == "pro")
    {
        # Retrieve the number of probes in promoters and genes
        NbrProbesPromoters <- read.csv(file.path(Origin, "matrices/M_Promoterss_NbrProbes.csv"),
                                   row.names = 1)
        print("Number of probes read")
        geneList <- merge(geneList, NbrProbesPromoters, by.x="Gene", by.y="row.names")
        #geneList$logFC <- -log10(geneList$adj.P.Val) * geneList$logFC
        geneList$logFC <- geneList$NbrProbe * geneList$logFC
        print(head(geneList))
    }
    
    # Create the dataframe needed for GSEA analysis
    corr <- bitr(geneList$Gene,
                 fromType  = "SYMBOL",
                 toType    = "ENTREZID",
                 OrgDb     = "org.Hs.eg.db")
    print(head(corr))
    d <- merge(corr, geneList, by.x = "SYMBOL", by.y = "Gene")[c("ENTREZID", "logFC")]
    print(head(d))
    geneList = d[,2]
    names(geneList) = as.character(d[,1])
    geneList = sort(geneList, decreasing = TRUE)
    print(head(geneList))

    #return(geneList)
    # Perform GSEA analysis
    try(kk2 <- gseKEGG(geneList       = geneList,
                       organism       = 'hsa',
                       keyType        = "ncbi-geneid",
                       minGSSize      = 2,
                       maxGSSize      = 500,
                       pvalueCutoff   = pvalueCutoff,
                       pAdjustMethod  = "BH",
                       verbose        = FALSE,
                       seed           = TRUE))

    # Return the dataframe
    if (exists("kk2"))
    {
        print(kk2)
        if (nrow(data.frame(kk2)) != 0)
        {
            Gse_Kegg <- setReadable(kk2,
                                    OrgDb   = org.Hs.eg.db,
                                    keyType ="ENTREZID")
            Gse_Kegg <- Gse_Kegg %>%
                arrange(desc(abs(NES)))

            slot(Gse_Kegg, "result")$NES_sign <- ifelse(slot(Gse_Kegg, "result")$NES < 0, 'hypomethylated', "hypermethylated")

            # Return data
            #print(map)
            ## Output images for all KEGG pathways found significantly enriched
            if (map)
            {                
                library("pathview")
                OutFilesF = file.path("~/ATM_Analysis/svn/analyse/results/KeggMaps/", comparison) 
                if(!dir.exists(OutFilesF))
                {
                    dir.create(OutFilesF, recursive = TRUE)
                }
                OutFilesF = file.path(OutFiles1, "Kegg/Gse/KeggMaps/", paste0(comparison,
                                                                              "_",
                                                                              feature)) 
                if(!dir.exists(OutFilesF))
                {
                    dir.create(OutFilesF, recursive = TRUE)
                }
                
                setwd(OutFilesF)
                print("folder created")
                get_kegg_plots <- function(GSEA, comparison)
                {
                    for (n in GSEA$ID)
                    {
                        pathview(gene.data = geneList,
                                 pathway.id = n,
                                 gene.idtype = "ENTREZID",#"SYMBOL",
                                 species = "hsa",
                                 kegg.dir = "~/ATM_Analysis/svn/analyse/results/KeggMaps/",
                                 limit = list(gene = max(abs(geneList)), cpd = 1))
                    }
                }
                get_kegg_plots(Gse_Kegg, comparison)
                print("function called")
            }

            if (!is.null(comparison))
            {
                slot(Gse_Kegg, "result")$comparison <- comparison
                return(Gse_Kegg)
            }
            else
            {
                return(Gse_Kegg)
            }

        } else {
            return(data.frame(matrix(ncol = 0, nrow = 0)))
        }
    }
}

    
#############   For Genes  ################
# Retrieve Kegg for Genes
A2.Genes.KeggGse <- getKeggGSE(data.genes.A2, comparison = "A2",
                               map = FALSE)
A0A2.Genes.KeggGse <- getKeggGSE(data.genes.A0A2, comparison = "A0A2",
                                 map = FALSE)
A2Inactive.Genes.KeggGse <- getKeggGSE(data.genes.A2Inactive, comparison = "A2Inactive",
                                       map = FALSE)


# Plot the comparison of the 4 analysis
GseKegg.all <- rbind(data.frame(A2.Genes.KeggGse),
                     data.frame(A0A2.Genes.KeggGse),
                     data.frame(A2Inactive.Genes.KeggGse))
print("binding done")
write.table(GseKegg.all, paste0(OutFiles, "KeggGse_Genes_according_comparisons.csv"), sep="\t")
print("written")


comment <- function(){
if (nrow(GseKegg.all) != 0){
    plot.genes.KeggGse <- dotplot_custom(GseKegg.all,
                                         main        = "Gene set enrichments analysis - Genes",
                                         comparison  = TRUE,
                                         orderAccordingTo = "signal")
    ggsaveDotPlot(filename = paste0(OutFiles, "KeggGse_comparisons_Genes.png"), plot = plot.genes.KeggGse)
}
}

#rm(A1.Genes.KeggGse, A2.Genes.KeggGse, A3.Genes.KeggGse, A4.Genes.KeggGse, all.Genes.KeggGse)
#rm(GseKegg.all, plot.genes.KeggGse, plot.genes.Kegg.all)





#############   For Promoters  ################

## For Promoters
# Retrieve Kegg for Promoters
A2.Pro.KeggGse <- getKeggGSE(data.pro.A2, comparison = "A2", map = FALSE)
A0A2.Pro.KeggGse <- getKeggGSE(data.pro.A0A2, comparison = "A0A2", map = FALSE)
A2Inactive.Pro.KeggGse <- getKeggGSE(data.pro.A2Inactive, comparison = "A2Inactive", map = FALSE)


print("all pro done")
# Plot the comparison of the 4 analysis
GseKegg.Pro <- rbind(data.frame(A2.Pro.KeggGse),
                     data.frame(A0A2.Pro.KeggGse),
                     data.frame(A2Inactive.Pro.KeggGse))
write.table(GseKegg.Pro, paste0(OutFiles, "KeggGse_Promoters_according_comparisons.csv"), sep=",")

print("and written")


comment <- function(){
if (nrow(GseKegg.Pro) != 0){
    plot.Pro.KeggGse <- dotplot_custom(GseKegg.Pro,
                                       main = "Gene set enrichments analysis - Promoters",
                                       comparison=TRUE,
                                       orderAccordingTo = "signal")
    ggsaveDotPlot(filename = paste0(OutFiles, "Comparisons_KeggGse_Pro.png"),
           plot = plot.Pro.KeggGse)
}
print("first dot plot done")

#rm(A1.Pro.KeggGse, A2.Pro.KeggGse, A3.Pro.KeggGse, A4.Pro.KeggGse, all.Pro.KeggGse)
#rm(GseKegg.Pro, plot.Pro.KeggGse, plot.Pro.Kegg.all)
}

save(A0A2.Pro.KeggGse,A2.Pro.KeggGse, A2Inactive.Pro.KeggGse,
     A0A2.Genes.KeggGse, A2.Genes.KeggGse, A2Inactive.Genes.KeggGse, 
     file = file.path(OutFiles, "AllKeggGseObjects.Rdata")
     )
print("object saved")


### Publication

if (nrow(A2.Pro.KeggGse) != 0){
    plot <- dotplot_custom_Paper(A2.Pro.KeggGse, orderAccordingTo = "GeneRatio", size = "NES", showSign=TRUE, showCategory = 50)
    ggsave(filename = paste0(OutFiles, "A2.Paper_KeggGse.png"),
           plot = plot,
           height = 45*0.25,
           width=12)
}

if (nrow(A0A2.Pro.KeggGse) != 0){
    plot <- dotplot_custom_Paper(A0A2.Pro.KeggGse, orderAccordingTo = "GeneRatio", size = "NES", showSign=TRUE, showCategory = 50)
    ggsave(filename = paste0(OutFiles, "A0A2.Paper_KeggGse.png"),
           plot = plot,
           height = 45*0.25,
           width=12)
}

if (nrow(A2Inactive.Pro.KeggGse) != 0){
    plot <- dotplot_custom_Paper(A2Inactive.Pro.KeggGse, orderAccordingTo = "GeneRatio", size = "NES", showSign=TRUE, showCategory = 50)
    ggsave(filename = paste0(OutFiles, "A2Inactive.Paper_KeggGse.png"),
           plot = plot,
           height = 45*0.25,
           width=12)
}




####### Test of lolliplots from GSEA results

plot <- lolliplot_custom_Paper(data = A2.Pro.KeggGse,
                               main              = NULL,
                               comparison        = FALSE,
                               orderAccordingTo  = "GeneRatio",
                               size              = "NES",
                               showSign          = T,
                               showCategory      = NULL,
                               sign              = NULL)

ggsave(file.path(OutFiles, "lolliplot_A2.png"), plot)

plot <- lolliplot_custom_Paper(data = A0A2.Pro.KeggGse,
                               main              = NULL,
                               comparison        = FALSE,
                               orderAccordingTo  = "GeneRatio",
                               size              = "NES",
                               showSign          = T,
                               showCategory      = NULL,
                               sign              = NULL)

ggsave(file.path(OutFiles, "lolliplot_A0A2.png"), plot)


plot <- lolliplot_custom_Paper(data = A2Inactive.Pro.KeggGse,
                               main              = NULL,
                               comparison        = FALSE,
                               orderAccordingTo  = "GeneRatio",
                               size              = "NES",
                               showSign          = T,
                               showCategory      = NULL,
                               sign              = NULL)

ggsave(file.path(OutFiles, "lolliplot_A2Inactive.png"), plot)



#grid.draw(KeggGse_comp)

print("all lolliplot custom")
KeggGse_comp <- lolliplot_custom_Paper(rbind(as.data.frame(A2.Pro.KeggGse),
                                             as.data.frame(A2Inactive.Pro.KeggGse),
                                             as.data.frame(A0A2.Pro.KeggGse)),
                                       comparison = TRUE,
                                       main              = NULL,
                                       orderAccordingTo  = "GeneRatio",
                                       size              = "NES",
                                       showSign          = T,
                                       showCategory      = NULL,
                                       sign              = NULL
                                       )

ggsave(file.path(OutFiles, "lolliplot_3_comparison.png"), KeggGse_comp, width = 17, height=20, units="in")


KeggGse_comp <- lolliplot_custom_Paper(rbind(as.data.frame(A2.Pro.KeggGse),
                                             as.data.frame(A2Inactive.Pro.KeggGse)),
                                       comparison = TRUE,
                                       main              = NULL,
                                       orderAccordingTo  = "GeneRatio",
                                       size              = "NES",
                                       showSign          = T,
                                       showCategory      = NULL,
                                       sign              = NULL
                                       )

ggsave(file.path(OutFiles, "lolliplot_2_comparison.png"), KeggGse_comp, width = 17, height=20, units="in")












## For Promoters
# Retrieve Kegg for Promoters
A2.Pro.KeggGse25 <- getKeggGSE(data.pro.A2, comparison = "A2", map = FALSE,
                               pvalueCutoff = 0.01, NbrProbes=TRUE)
print("A0 done")
A0A2.Pro.KeggGse25 <- getKeggGSE(data.pro.A0A2, comparison = "A0A2", map = FALSE,
                                 pvalueCutoff = 0.01, NbrProbes=TRUE)
print("A0A2 done")
A2Inactive.Pro.KeggGse25 <- getKeggGSE(data.pro.A2Inactive, comparison = "A2Inactive", map = FALSE,
                                 pvalueCutoff = 0.01, NbrProbes=TRUE)
print("A2Inactive done")

print("all pro done")
# Plot the comparison of the 4 analysis
GseKegg.Pro <- rbind(data.frame(A2.Pro.KeggGse25),
                     data.frame(A0A2.Pro.KeggGse25),
                     data.frame(A2Inactive.Pro.KeggGse25))
write.table(GseKegg.Pro, paste0(OutFiles, "Nbr_KeggGse0_01_Promoters_according_comparisons.csv"), sep=",")
print("and written")

if (nrow(GseKegg.Pro) != 0){
    plot.Pro.KeggGse <- dotplot_custom(GseKegg.Pro,
                                       main = "Gene set enrichments analysis - Promoters",
                                       comparison=TRUE,
                                       orderAccordingTo = "signal")
    ggsaveDotPlot(filename = paste0(OutFiles, "Nbr_Comparisons_KeggGse0_01_Pro.png"),
           plot = plot.Pro.KeggGse)
}
print("first dot plot done")

#rm(A1.Pro.KeggGse, A2.Pro.KeggGse, A3.Pro.KeggGse, A4.Pro.KeggGse, all.Pro.KeggGse)
#rm(GseKegg.Pro, plot.Pro.KeggGse, plot.Pro.Kegg.all)


save(A0A2.Pro.KeggGse25, A2Inactive.Pro.KeggGse25, A2.Pro.KeggGse25,
     file = file.path(OutFiles, "Nbr_AllKeggGseObjects0_01.Rdata")
     )
print("object saved")


### Publication
comment <- function(){
if (nrow(A2.Pro.KeggGse25) != 0){
    plot <- dotplot_custom_Paper(A2.Pro.KeggGse25, orderAccordingTo = "GeneRatio", size = "NES", showSign=TRUE, showCategory = 50)
    ggsave(filename = paste0(OutFiles, "A2.Paper_KeggGse0_01.png"),
           plot = plot,
           height = 45*0.25,
           width=12)
}

if (nrow(A0A2.Pro.KeggGse25) != 0){
    plot <- dotplot_custom_Paper(A0A2.Pro.KeggGse25, orderAccordingTo = "GeneRatio", size = "NES", showSign=TRUE, showCategory = 50)
    ggsave(filename = paste0(OutFiles, "A0A2.Paper_KeggGse0_01.png"),
           plot = plot,
           height = 45*0.25,
           width=12)
}

if (nrow(A2Inactive.Pro.KeggGse25) != 0){
    plot <- dotplot_custom_Paper(A2Inactive.Pro.KeggGse25, orderAccordingTo = "GeneRatio", size = "NES", showSign=TRUE, showCategory = 50)
    ggsave(filename = paste0(OutFiles, "A2Inactive.Paper_KeggGse0_01.png"),
           plot = plot,
           height = 45*0.25,
           width=12)
}
}



KeggGse_comp <- lolliplot_custom_Paper(rbind(as.data.frame(A2.Pro.KeggGse25),
                                             as.data.frame(A2Inactive.Pro.KeggGse25)),
                                             #as.data.frame(AVUS.Pro.KeggGse),
                                             #as.data.frame(A2Inactive.Pro.KeggGse)),
#                                       categories = label,
                                       comparison = TRUE,
                                       main              = NULL,
                                       orderAccordingTo  = "GeneRatio",
                                       size              = "NES",
                                       showSign          = T,
                                       showCategory      = NULL,
                                       sign              = NULL
                                       )

ggsave(file.path(OutFiles, "Nbr_lolliplot_2_comparison0_01.png"), KeggGse_comp, width = 17, height=20, units="in")










comment <- function(){
# By removing GPR161

## For Promoters
# Retrieve Kegg for Promoters
A2.Pro.KeggGseGPR161 <- getKeggGSE(data.pro.A2[data.pro.A2$Gene != "GPR161", ],
                                   comparison = "A2", map = FALSE)
print("A0 done")
A0A2.Pro.KeggGseGPR161 <- getKeggGSE(data.pro.A2[data.pro.A2$Gene != "GPR161", ],
                                 comparison = "A0A2", map = FALSE)
print("A0A2 done")
A2Inactive.Pro.KeggGseGPR161 <- getKeggGSE(data.pro.A2[data.pro.A2$Gene != "GPR161", ],
                                       comparison = "A2Inactive", map = FALSE,)
print("A2Inactive done")

print("all pro done")
# Plot the comparison of the 4 analysis
GseKegg.Pro <- rbind(data.frame(A2.Pro.KeggGseGPR161),
                     data.frame(A0A2.Pro.KeggGseGPR161),
                     data.frame(A2Inactive.Pro.KeggGseGPR161))
write.table(GseKegg.Pro, paste0(OutFiles, "KeggGseGPR161_Promoters_according_comparisons.csv"), sep=",")
print("and written")

}




#############   For Promoters  ################

## For Promoters
# Retrieve Kegg for Promoters
A2.Pro.KeggGseNbr <- getKeggGSE(data.pro.A2, comparison = "A2", map = FALSE,
                                NbrProbes = TRUE)
print("A0 done")
A0A2.Pro.KeggGseNbr <- getKeggGSE(data.pro.A0A2, comparison = "A0A2", map = FALSE,
                                  NbrProbes = TRUE)
print("A0A2 done")
A2Inactive.Pro.KeggGseNbr <- getKeggGSE(data.pro.A2Inactive, comparison = "A2Inactive", map = FALSE,
                                    NbrProbes = TRUE)
print("A2Inactive done")


# Plot the comparison of the 4 analysis
GseKegg.Pro <- rbind(data.frame(A2.Pro.KeggGseNbr),
                     data.frame(A0A2.Pro.KeggGseNbr),
                     data.frame(A2Inactive.Pro.KeggGseNbr))
write.table(GseKegg.Pro, paste0(OutFiles, "Nbr_KeggGse_Promoters_according_comparisons.csv"), sep=",")
print("and written")



KeggGse_comp <- lolliplot_custom_Paper(rbind(as.data.frame(A2.Pro.KeggGseNbr),
                                             as.data.frame(A2Inactive.Pro.KeggGseNbr)),
                                       comparison = TRUE,
                                       main              = NULL,
                                       orderAccordingTo  = "GeneRatio",
                                       size              = "NES",
                                       showSign          = T,
                                       showCategory      = NULL,
                                       sign              = NULL,
                                       ySize             = 10,
                                       )

ggsave(file.path(OutFiles, "Nbr_KeggGse_lolliplot_2_comparison.png"), KeggGse_comp, width = 10, height=14, units="in")




stop()


load(file.path(Origin, "/matrices/M_Funnorm_Normalisation.RData"))
load(file.path(Origin, "/matrices/Mapping.RData"))





preComputePheatmap <- function(data,
                               annotations,
                               colAnno = c("Supposed biallelic inactivation", "Sample_group"),
                               method_dist,
                               method_clust,
                               cutree_cols = 2,
                               cutree_rows = 2,
                               showRownames = F,
                               showColnames = F,
                               fontsize_row = 12)
{
    color <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100)
    
    hc.features <- hclust(as.dist(1-cor(t(data), method=method_dist)), method=method_clust)
    hc.samples <- hclust(as.dist(1-cor(data, method=method_dist)), method=method_clust)

    print(table(annotations[['label']]))
    annotations <- annotations %>% mutate(Inactive = recode(Inactive,
                                                            "active ATM" = "No",
                                                            "Control" = "No",
                                                            "inactive ATM" = "Yes"),
                                          label = recode(ClinVar,
                                                         "Likely pathogenic" = "PV",
                                                         "Pathogenic" = "PV",
                                                         "Pathogenic/Likely pathogenic" = "PV",
                                                         "Uncertain significance" = "VUS")
)
    
    annotations[["Supposed biallelic inactivation"]] <- annotations[["Inactive"]]
    annotations$"Supposed biallelic inactivation" <- factor(annotations$"Supposed biallelic inactivation", levels=c("Yes", "No"))

    # Sample_Group
    annotations[['Sample_Group']] <- paste0(annotations[['Sample_Group']], '_', annotations[['label']])

    annotations <- annotations %>% mutate(Sample_Group = recode(Sample_Group,
                                                                "Non.ATM_" = "Non-ATM",
                                                                "ATM_PV"="ATM PV",
                                                                "ATM_VUS"="ATM VUS"))
    # Grade
    annotations["Grade"][annotations["Grade"] == ""] <- NA
    annotations <- annotations %>% dplyr::mutate(Grade = replace_na(Grade, "Unknown"))

    # Subtype
    print(table(annotations$subtype))
    annotations["subtype"][annotations["subtype"] == ""] <- NA
    annotations["subtype"][annotations["subtype"] == "Missing"] <- NA
        annotations["subtype"][annotations["subtype"] == "No record"] <- NA
    annotations <- annotations %>% dplyr::mutate(subtype = replace_na(subtype, "Unknown"))
#    annotations <- annotations %>% mutate("Sample_Group" = recode(Sample_Group,
#                                                                  "ATM PV"="ATM_PV",
#                                                                  "ATM VUS"="ATM_VUS"))
    print(table(annotations$subtype))
    annoCol <- list(
        "Sample_Group" = c("ATM PV"="#00DAE0", "ATM VUS"="#FDD513",
                           "Non-ATM"="#E199FF"),
        "Supposed biallelic inactivation" = c(Yes="#FF9289", No="#96CA00"),
        "Grade" = c("I"="yellow", "II"="orange", "III"="red", "Unknown"="white"),
        "subtype" = c("Luminal"="light blue", "Luminal A"="blue", "Luminal B"="dark blue",
                      "Luminal B/HER2+"="purple", "Triple negative"="black", "HER2+"="green",
                      "Unknown"="white")
    )

    newnames <- lapply(
          rownames(data),
          function(x) bquote(italic(.(x))))    
            
    heatmap <- pheatmap::pheatmap(data,
                                  treeheight_row     = 0,
                                  annotation_col     = annotations[colAnno],
                                  show_colnames      = showColnames,
                                  show_rownames      = showRownames,
                                  annotation_colors =  annoCol,
                       #custering_distance_cols= "correlation",
                       #custering_distance_rows= "correlation",
                                  scale              = "row",
                                  cluster_rows       = hc.features,
                                  cluster_cols       = hc.samples,
                                  cutree_cols        = cutree_cols,
                                  cutree_rows        = cutree_rows,
                                  color              = color,
                                  labels_row = as.expression(newnames),
                                  fontsize_row = fontsize_row)
    return(heatmap)
}



# Vizualise heatmap for DM genes
OutFilesGenesAnalysis = file.path(Origin, "genesAnalysis")

# all promoters DF in A2
H1 <- preComputePheatmap(data=data_Promoters[pro.A2$Gene,
                                      annotations_A2$Sample_Name],
                   annotations = annotations_A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   colAnno = c("Grade","subtype", "Supposed biallelic inactivation", "Sample_Group"),
                   )

# all promoters DF in A2Inactive
H2 <- preComputePheatmap(data=data_Promoters[pro.A2Inactive$Gene,
                                      annotations_A2Inactive$Sample_Name],
                   annotations = annotations_A2Inactive,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2
                   )

# all promoters DF in A0A2
H3 <- preComputePheatmap(data=data_Promoters[pro.A0A2$Gene,
                                      annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2
                   )

library(patchwork)
all <- ggplotify::as.ggplot(H1) + ggplotify::as.ggplot(H2) + ggplotify::as.ggplot(H3) +
    plot_layout(ncol=2) + plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(OutFilesGenesAnalysis, "Heatmaps_DM_promoters.png"),
       plot     = all,
       width = 22.4,
       height = 11.725,
       units = "in")




# all promoters DF in A2 with all samples
H1.2 <- preComputePheatmap(data=data_Promoters[pro.A2$Gene,
                                      annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2
                   )

library(patchwork)
all <- ggplotify::as.ggplot(H1) + ggplotify::as.ggplot(H1.2) + 
    plot_layout(ncol=2) + plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(OutFilesGenesAnalysis, "Heatmaps_DM_A2_with_VUS.png"), 
       plot     = all,
       width = 22.4,
       height = 11.725,
       units = "in")







# all promoters DF in A2 with all samples
H2.2 <- preComputePheatmap(data=data_Promoters[pro.A2Inactive$Gene,
                                      annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2
                   )

library(patchwork)
all <- ggplotify::as.ggplot(H1) + ggplotify::as.ggplot(H1.2) + 
    plot_layout(ncol=2) + plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(OutFilesGenesAnalysis, "Heatmaps_DM_A2_with_VUS.png"), 
       plot     = all,
       width = 22.4,
       height = 11.725,
       units = "in")







genesCommon <- intersect(pro.A0A2$Gene, intersect(pro.A2$Gene, pro.A2Inactive$Gene))


preComputePheatmap(data=data_Promoters[,
                                       annotations_A0A2[annotations_A0A2$Sample_Group == "ATM",]$Sample_Name],
                   annotations = annotations_A0A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2
                   )










heat <- t(scale(t(data_Promoters[pro.A2Inactive$Gene,
                                 annotations_A2Inactive$Sample_Name])))
annotations_A2Inactive$"ATM activity" <- annotations_A2Inactive$Inactive
pheatmap::pheatmap(heat, treeheight_row = 0,
                   annotation_col = annotations_A2Inactive[c("ATM activity", "Sample_Group")],
                   scale = "row",
                   color= color,
                   custering_distance = "pearson",
                   clustering_method = "ward.D2",
                   cutree_cols = 2,
                   show_rownames=F,
                   show_colnames=F)






# Genes selected with 






umap_func(data_Promoters[pro.A2$Gene,
                                      annotations_A2$Sample_Name],
          sampleTable = annotations_A0A2,
          highlightedVarName = "Sample_Group",
          highlightedVar = "Sample_Group",
          metric = "pearson2",          
          seed=0)




annotations_A0A2$ER_status <- as.character(annotations_A0A2$ER_status)

annotations_A0A2$ER_status2 <- annotations_A0A2$ER_status
annotations_A0A2[is.na(annotations_A0A2$ER_status2), ]$ER_status2 <- "NA"
annotations_A0A2[annotations_A0A2$ER_status2 == "1", ]$ER_status2 <- "Positive"
annotations_A0A2[annotations_A0A2$ER_status2 == "0", ]$ER_status2 <- "Negative"
annotations_A0A2[annotations_A0A2$ER_status2 == "NA", ]$ER_status2 <- NA


umap <- umap_func_publiv2(matrix                      = data_Promoters[pro.A0A2$Gene,
                                      annotations_A0A2$Sample_Name],
                          sampleTable                 = annotations_A0A2,
                          highlightedVar              = "ER_status2",
                          highlightedVarName          = "ER status",
                          shape                       = "Inactive",
                          colors                      = colors,
                          fontsize                    = 14,
                          pointSize                   = 2,
                          metric                      = "pearson2",
                          title                       = NULL,#"UMAP using the 10.000 most variable probes",
                          colCond = "Sample_Group",
                          valCond = "ATM",
                          seed = 0) 




# Visualise heatmap for genes selected by ML methods

pathways <- read.csv(file.path(Origin, "Interpretation/Kegg/Gse/Nbr_KeggGse_Promoters_according_comparisons.csv"))
A2.path <- pathways[pathways$comparison == "A2", ]
A2.path <- ParseEnrichResults(A2.path, column = "core_enrichment")
pathways <- c("Amyotrophic lateral sclerosis", "Hepatocellular carcinoma",
              "Huntington disease", "Prion disease", "Ribosome",
              "Shigellosis")

genes <- c(A2.path[A2.path$Description == pathways[1], ]$core_enrichment,
           A2.path[A2.path$Description == pathways[2], ]$core_enrichment,
           A2.path[A2.path$Description == pathways[3], ]$core_enrichment,
           A2.path[A2.path$Description == pathways[4], ]$core_enrichment,
           A2.path[A2.path$Description == pathways[5], ]$core_enrichment)
genes <- genes[!duplicated(genes)]


preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A2$Sample_Name],
                   annotations = annotations_A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2
                   )



color <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)    
heat <- t(scale(t(data_Promoters[genes,
                       annotations_A2$Sample_Name])))
pheatmap::pheatmap(heat,
                   treeheight_row      = 0,
                   annotation_col      = annotations_A2["Sample_Group"],
                   custering_distance_cols= "correlation",
                   custering_distance_rows= "correlation",
                   clustering_method  = "ward.D2",
                   scale               = "row",
                   cutree_cols         = 2,
                   color               = color,
                   width               = 8,
                   height              = 6)



genes <- c('API5P2','CAMKK2','CSTB','EIF3B','FABP6','GNPTAB','H3-3A','HK2','IMPDH1P10','INTS6P1','LGALS8-AS1','LGALSL','LINC01600','MFHAS1','MIR6892','PABIR1','PKN2-AS1','PTDSS2','PTPRVP','RNU6-72P','RPL17','RPL17-C18orf32','SNORD73A','SPG21','ST8SIA6-AS1','SUSD3','TRMT112P6','UBE2R2','ZER1')

umap <- umap_func_publiv2(matrix                      = data_Promoters[genes,
                                                                       annotations_A2$Sample_Name],
                          sampleTable                 = annotations_A2,
                          highlightedVar              = "Sample_Group",
                          highlightedVarName          = "Sample_Group",
                          shape                       = "Inactive",
                          colors                      = brewer.pal(n = 3, name = "Dark2"),
                          fontsize                    = 14,
                          pointSize                   = 3,
                          metric                      = "pearson2",
                          title                       = NULL)

ggsave("~/ATM_Analysis/svn/analyse/result/FinalMethylation/test.png", test)


genes <- c('API5P2','CAMKK2','CSTB','EIF3B','FABP6','GNPTAB','H3-3A','HK2','IMPDH1P10','INTS6P1','LGALS8-AS1','LGALSL','LINC01600','MFHAS1','MIR6892','PABIR1','PKN2-AS1','PTDSS2','PTPRVP','RN7SKP114','RNU6-72P','RPL17','RPL17-C18orf32','SNORD73A','SPG21','ST8SIA6-AS1','SUSD3','TRMT112P6','UBE2R2','ZER1')

genes <- c('API5P2','CAMKK2','CSTB','EIF3B','FABP6','GNPTAB','H3-3A','HK2','IMPDH1P10','INTS6P1','LGALS8-AS1','LGALSL','LINC01600','MFHAS1','MIR6892','PABIR1','PKN2-AS1','PTDSS2','PTPRVP','RN7SKP114','RNU6-72P','RPL17','RPL17-C18orf32','SNORD73A','SPG21','ST8SIA6-AS1','SUSD3','TRMT112P6','UBE2R2','ZER1')


umap <- umap_func_publiv2(matrix                      = data_Promoters[genes,
                                                                       annotations_A2$Sample_Name],
                          sampleTable                 = annotations_A2,
                          highlightedVar              = "Sample_Group",
                          highlightedVarName          = "Sample_Group",
                          shape                       = "Inactive",
                          colors                      = brewer.pal(n = 3, name = "Dark2"),
                          fontsize                    = 14,
                          pointSize                   = 3,
                          metric                      = "pearson2",
                          title                       = NULL)

ggsave("~/ATM_Analysis/svn/analyse/result/FinalMethylation/test.png", test)




genes <- c("SNORD91B", "SCAPER")

genes <- c('IGKV1OR2-108','INTS6P1','MALT1-AS1','PDIA3P2','PKN2-AS1','POLR2L','PRH1','PRR4','RASGEF1C','RPL21P76','RPL36AP30','SCAPER','SLC7A5P2','SMG1P3','SMIM10L1','SNORD91A','SNORD91B')

genes <- c('AMPD3','ARF4','CSTB','HNRNPA1P49','HSPE1P10','IGKV1OR2-108','INTS6P1','LINC00327','MALT1-AS1','MIR1266','MIR6892','PDIA3P2','PKN2-AS1','POLR2L','PRH1','PRR4','PTDSS2','RASGEF1C','RPL21P76','RPL36AP30','SCAPER','SLC7A5P2','SMG1P3','SMIM10L1','SNORA14A','SNORD88C','SNORD91A','SNORD91B','TAS2R14')


umap <- umap_func_publiv2(matrix                      = data_Promoters[genes,
                                                                       annotations_A2$Sample_Name],
                          sampleTable                 = annotations_A2,
                          highlightedVar              = "Sample_Group",
                          highlightedVarName          = "Sample_Group",
                          shape                       = "Inactive",
                          colors                      = brewer.pal(n = 3, name = "Dark2"),
                          fontsize                    = 14,
                          pointSize                   = 3,
                          metric                      = "pearson2",
                          title                       = NULL)

ggsave("~/ATM_Analysis/svn/analyse/result/FinalMethylation/test.png", test)




genes <- c("EIF3B", "INTS6P1", "PKN2-AS1")

preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   colAnno = c("ATM activity", "Sample_Group"),
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames = TRUE
                   )






umap_func_publiv2(matrix       = data_Promoters[genes,
                                                annotations_A0A2$Sample_Name],
                  sampleTable                 = annotations_A0A2,
                  highlightedVar              = "Sample_Group",
                  highlightedVarName          = "Sample_Group",
                  shape                       = "Inactive",
                  colors                      = brewer.pal(n = 3, name = "Dark2"),
                  fontsize                    = 14,
                  pointSize                   = 3,
                  metric                      = "pearson2",
                  title                       = NULL,
                  n_neighbors = 5)





# All pmoters after feature selection
preComputePheatmap(data=data_Promoters[unname(unlist(all.features)),
                                       annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   colAnno = c("ATM activity", "Sample_Group"),
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames = TRUE
                   )


# All promoters after feature selection for each sensitivity analysis

genes <- unname(unlist(all.features[names(all.features)[grepl("^A2_", names(all.features))]]))
genes <- genes[!duplicated(genes)]
H1 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A2$Sample_Name],
                   annotations = annotations_A2,
                   colAnno = c("ATM activity", "Sample_Group"),
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames = TRUE,
                   fontsize_row = 7.5)


genes <- unname(unlist(all.features[names(all.features)[grepl("^A0A2_", names(all.features))]]))
genes <- genes[!duplicated(genes)]
H3 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   colAnno = c("ATM activity", "Sample_Group"),
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames = TRUE,
                   fontsize_row = 8)

genes <- unname(unlist(all.features[names(all.features)[grepl("^A2Inactive_", names(all.features))]]))
genes <- genes[!duplicated(genes)]
H2 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A2Inactive$Sample_Name],
                   annotations = annotations_A2Inactive,
                   colAnno = c("ATM activity", "Sample_Group"),
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames = TRUE,
                   fontsize_row = 8)



library(patchwork)
all <- ggplotify::as.ggplot(H1) + ggplotify::as.ggplot(H2) + ggplotify::as.ggplot(H3) +
    plot_layout(ncol=2) + plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(Origin, "Heatmaps_FS_promoters.png"),
       plot     = all,
       width = 22.4,
       height = 11.725,
       units = "in")



# Stable genes found at least twice for each sensitivity analysis

genes <- c("EIF3B", "INTS6P1", "PKN2-AS1", "PDIA3P2", "RASGEF1C", "LGALS8-AS1")
P1 <- preComputePheatmap(data=data_Promoters[genes,
                                             annotations_A2$Sample_Name],
                   annotations = annotations_A2,
                   colAnno = c("ATM activity", "Sample_Group"),
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames = TRUE
                   )

P3 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   colAnno = c("ATM activity", "Sample_Group"),
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames = TRUE
                   )

genes <- c("EIF3B", "INTS6P1", "PKN2-AS1", "CSTB", "SNORA14A")
P2 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A2Inactive$Sample_Name],
                   annotations = annotations_A2Inactive,
                   colAnno = c("ATM activity", "Sample_Group"),
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames = TRUE
                   )



library(patchwork)
all <- ggplotify::as.ggplot(P1) + ggplotify::as.ggplot(P2) + ggplotify::as.ggplot(P3) +
    plot_layout(ncol=3) + plot_annotation(tag_levels = 'a')

ggsave(filename = file.path(Origin, "Heatmaps_stable-FS_promoters.png"),
       plot     = all,
       width = 22.4,
       height = 8,#11.725,
       units = "in")





multipleUMAP <- function(all.features)
{

    listPlots <- list()
    for (i in seq_along(all.features))
    {
        print(names(all.features[i]))
        n <- sub("^(.*?)_.*", "\\1", names(all.features[i]))
        print(n)
        anno = eval(parse(text=paste("annotations_", n, sep="")))
        
        p <- umap_func_publiv2(matrix       = data_Promoters[
                                   all.features[[names(all.features)[i]]],
                                   anno$Sample_Name],
                               sampleTable                 = anno,
                               highlightedVar              = "Sample_Group",
                               highlightedVarName          = "Sample_Group",
                               shape                       = "Inactive",
                               colors                      = brewer.pal(n = 3, name = "Dark2"),
                               fontsize                    = 14,
                               pointSize                   = 3,
                               metric                      = "pearson2",
                               title                       = names(all.features[i]),
                               n_neighbors = 5)
        
        listPlots[[i]] <- p
    }
    return(listPlots)
}


plot_list <- multipleUMAP(all.features)

p1 <- patchwork::wrap_plots(plot_list,
                            nrow = 3, ncol = 3)

p1 <- p1 + theme(legend.position = "right") + patchwork::plot_layout(guides = "collect")

ggsave("~/ATM_Analysis/svn/analyse/results/FinalMethylation/feature_selection_umapFinal.png", p1, width=25.6, height=13.4, units="in")



unname(unlist(all.features))[unname(unlist(all.features)) %in% oncoKB]



d <- umap_func_publiv3(matrix       = data_Promoters[
                           all.features[[names(all.features)[i]]],
                           annotations_A2$Sample_Name],
                       sampleTable                 = annotations_A2,
                       highlightedVar              = "Sample_Group",
                       highlightedVarName          = "Sample_Group",
                       shape                       = "Inactive",
                       colors                      = brewer.pal(n = 3, name = "Dark2"),
                       fontsize                    = 14,
                       pointSize                   = 3,
                       metric                      = "pearson2",
                       title                       = names(all.features[i]),
                       n_neighbors = 5,
                       returnData = T)




p <- plot_ly(x=d$UMAP1, y=d$UMAP2, z=d$UMAP3, type="scatter3d", mode="markers", color=d$Sample_Group)

htmlwidgets::saveWidget(as_widget(p), "~/ATM_Analysis/svn/analyse/results/FinalMethylation/3D_umap.html" )



#############################
## After feature selection ##
#############################

library("readxl")


readLists <- function(path)
{
    headers <- c("Logistic Regression", "Random Forest", "XGBoost")
    lines <- readLines(path)
    print(lines)
#    rm(dict, subList, listString, name)
    dict <- list()
    for (i in 1:length(lines)){
        if (i == length(lines)){
            print("I am at the end")
            dict[[name]] <- unlist(subList)                    
            subList <- list()
        } else if (lines[i] != ""){
            if (lines[i] %in% headers) {
                if (exists("name")){
                    dict[[name]] <- unlist(subList)                    
                    subList <- list()
                } else {
                    print("Name do not exists yet")
                    name <- lines[i]
                subList <- list()
                }
                name <- lines[i]
            } else {
                # Supprimer les crochets et les guillemets simples
                listString <- gsub("\\[|\\]|'", "", lines[i])
                # Diviser la chaîne de caractères en une liste
                listString <- strsplit(listString, ", ")[[1]]
                subList <- c(subList, listString)
            }
        }
    }
    return(dict)
}


path <- file.path(Origin, "Classifier/final_analyses/A2_All_features_selected.txt")
A2lists <- readLists(path)



path <- file.path(Origin, "Classifier/final_analyses/A2_other_tests_validation_3Final_features_selected.xlsx")
A2lists <- read_excel(path)





fill_merged <- function(dat, columns.as.vector){
    require('readxl')
    ## Check if column names are provided as strings
    if(!is.character(columns.as.vector)){
        stop("Column names must be provided as string or vector of strings of class character")
    }
    ## Go through the columns
    for(column in columns.as.vector){
        ## Check if the column name matches with dat column names
        if (!column %in% names(dat)){
            stop(paste0('Column <', column, '> cannot be found in the data frame'))
        }
        ## Get value of each row
        for(n in 1:nrow(dat)){
            ## Check if it is empty
            if(dat[[column]][n] == '' || is.na(dat[[column]][n])){
                ## If it is the row 1, stop with Error
                if(n == 1){
                    stop(paste0("Row 1 of column <", column,
                                    "> has empty values. Check your data."))
                }
                else{
                    dat[[column]][n] <- dat[[column]][n - 1]
                }
            }
        }
    }
    return(dat)
}

library(ggplot2)  # Load ggplot2 for violin plots

features.A2 <- read_xlsx(file.path(Origin, '/Classifier/A2_features_selected.xlsx'))
colnames(features.A2) <- paste0("A2_", colnames(features.A2), sep="") 
features.A2Inactive <- read_xlsx(file.path(Origin, '/Classifier/A2Inactive_features_selected.xlsx'))
colnames(features.A2Inactive) <- paste0("A2Inactive_", colnames(features.A2Inactive), sep="")

features.A0A2 <- read_xlsx(file.path(Origin, '/Classifier/A0A2_features_selected.xlsx'))
colnames(features.A0A2) <- paste0("A0A2_", colnames(features.A0A2), sep="")


all.features <- plyr::rbind.fill(features.A2, features.A2Inactive, features.A0A2)
all.features <- as.list(all.features)
all.features <- lapply(all.features, function(x) x[!is.na(x)])



library(UpSetR)

upset <- ggplotify::as.ggplot(upset(fromList(all.features),# sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"),
      sets = rev(c("A2_Logistic Regression", "A2Inactive_Logistic Regression", "A0A2_Logistic Regression",
                   "A2_Random Forest", "A2Inactive_Random Forest", "A0A2_Random Forest",
                   "A2_XGBoost", "A2Inactive_XGBoost", "A0A2_XGBoost")),
      keep.order = T, 
      empty.intersections = NULL,
      nsets = 6,
      sets.bar.color = "#56B4E9",
      order.by = "freq",
      text.scale = 2)
)

ggsave(filename = file.path(Origin, "/Classifier/upset.png"), upset,
       width = 50,
       height = 30,
       units = "cm",
       dpi=300)
                         



from_list <- function(dataframe) {
    members = unique(unlist(as.list(dataframe), use.names=F))
    members = members[!is.na(members)] 
    return(data.frame(sapply(dataframe, function(set) members %in% set),
                      row.names = members))
}


matrix = from_list(all.features)
matrix$gene_name = rownames(matrix)
head(matrix)


library(ComplexUpset)
upset(
    matrix,
    intersect=colnames(matrix)[1:length(colnames(matrix)) - 1],
    base_annotations=list(
        'Intersection size'=(
            intersection_size(
                bar_number_threshold=1,
                color='grey9',
                fill='grey80'
            )
            + geom_text(
                  mapping=aes(label=gene_name),
                  size = 2.8,
                  position=position_stack(),
                  na.rm=TRUE,
                  vjust=2.5
              )
        )
    ),
    width_ratio=0.15,
    height_ratio=1/4
)



classif <- read_xlsx(file.path(Origin, '/Classifier/pathways/A2_Pathways_evaluation_feature_selection.xlsx'))

classif <- fill_merged(classif, "FS method")
classif <- fill_merged(classif, "Classifier")
classif <- classif[c("FS method", "Classifier", "test_Matthews")]
colnames(classif)[colnames(classif) == "FS method"] <- "pathway"
classif <- classif %>%
    group_by(pathway, Classifier) %>%
    summarise(mean = mean(test_Matthews),
              std = sd(test_Matthews))

plot <- ggplot(
    classif[classif$Classifier %in% c("Random Forest",
                                      "XGBoost", "Support Vector Machine",
                                      "Logistic Regression"), ],
    aes(x = pathway, y = mean)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    geom_errorbar(aes(x = pathway,
                      ymin = mean - std,
                      ymax = mean + std)) +
    geom_hline(yintercept = 0.80, linetype=2, color = "red") + 
    ggtitle("Distribution of Matthews scores across Cross-Validation Folds") +
    xlab("Pathway") + ylab("Matthews score") + theme_bw() +
    theme(text=element_text(size=20),
          #axis.text.x = element_text(size=15),#angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) + 
    facet_wrap(~Classifier, as.table = FALSE, nrow=1) + coord_flip()

ggsave(filename = file.path(Origin, "/Classifier/pathways/A2_MCCS_pathways.png"), plot,
       width = 60,
       height = 50,
       units = "cm",
       dpi=300)


classif <- read_xlsx(file.path(Origin, '/Classifier/pathways/A2Inactive_Pathways_evaluation_feature_selection.xlsx'))

classif <- fill_merged(classif, "FS method")
classif <- fill_merged(classif, "Classifier")
classif <- classif[c("FS method", "Classifier", "test_Matthews")]
colnames(classif)[colnames(classif) == "FS method"] <- "pathway"
classif <- classif %>%
    group_by(pathway, Classifier) %>%
    summarise(mean = mean(test_Matthews),
              std = sd(test_Matthews))

plot <- ggplot(
    classif[classif$Classifier %in% c("Random Forest",
                                      "XGBoost", "Support Vector Machine",
                                      "Logistic Regression"), ],
    aes(x = pathway, y = mean)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    geom_errorbar(aes(x = pathway,
                      ymin = mean - std,
                      ymax = mean + std)) +
    geom_hline(yintercept = 0.80, linetype=2, color = "red") + 
    ggtitle("Distribution of Matthews scores across Cross-Validation Folds") +
    xlab("Pathway") + ylab("Matthews score") + theme_bw() +
    theme(text=element_text(size=20),
          #axis.text.x = element_text(size=15),#angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) + 
    facet_wrap(~Classifier, as.table = FALSE, nrow=1) + coord_flip()

ggsave(filename = file.path(Origin, "/Classifier/pathways/A2Inactive_MCCS_pathways.png"), plot,
       width = 60,
       height = 50,
       units = "cm",
       dpi=300)



# All promoters
classif <- read_xlsx(file.path(Origin, '/Classifier/A2_evaluation_several_ML.xlsx'))

classif <- fill_merged(classif, "Classifier")
classif <- classif[c("Classifier", "Fold", "test_Matthews")]


classif <- classif %>%
    group_by(Classifier) %>%
    summarise(mean = mean(test_Matthews),
              std = sd(test_Matthews))

head(classif)

plot <- ggplot(
    classif[classif$Classifier %in% c("Random Forest",
                                      "XGBoost", "Support Vector Machine",
                                      "Logistic Regression"), ],
    aes(x = Classifier, y = mean, fill = Classifier)) +
    geom_bar(stat = "identity", ) +
    geom_errorbar(aes(x = Classifier,
                      ymin = mean - std,
                      ymax = mean + std)) +
    geom_hline(yintercept = 0.80, linetype=2, color = "red") + 
    ggtitle("Distribution of Matthews scores across Cross-Validation Folds") +
    xlab("Classifier") + ylab("Matthews score") + theme_bw() +
    theme(text=element_text(size=20),
          #axis.text.x = element_text(size=15),#angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) #+ 
#    facet_wrap(~Classifier, as.table = FALSE, nrow=1) + coord_flip()

ggsave(filename = file.path(Origin, "/Classifier/A2_MCCS_genome-wide.png"), plot,
       width = 60,
       height = 50,
       units = "cm",
       dpi=300)







classif <- read_xlsx(file.path(Origin, '/Classifier/A2Inactive_evaluation_several_ML.xlsx'))

classif <- fill_merged(classif, "Classifier")
classif <- classif[c("Classifier", "Fold", "test_Matthews")]


classif <- classif %>%
    group_by(Classifier) %>%
    summarise(mean = mean(test_Matthews),
              std = sd(test_Matthews))

head(classif)

plot <- ggplot(
    classif[classif$Classifier %in% c("Random Forest",
                                      "XGBoost", "Support Vector Machine",
                                      "Logistic Regression"), ],
    aes(x = Classifier, y = mean, fill = Classifier)) +
    geom_bar(stat = "identity", ) +
    geom_errorbar(aes(x = Classifier,
                      ymin = mean - std,
                      ymax = mean + std)) +
    geom_hline(yintercept = 0.80, linetype=2, color = "red") + 
    ggtitle("Distribution of Matthews scores across Cross-Validation Folds") +
    xlab("Classifier") + ylab("Matthews score") + theme_bw() +
    theme(text=element_text(size=20),
          #axis.text.x = element_text(size=15),#angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) #+ 
#    facet_wrap(~Classifier, as.table = FALSE, nrow=1) + coord_flip()

ggsave(filename = file.path(Origin, "/Classifier/A2Inactive_MCCS_genome-wide.png"), plot,
       width = 60,
       height = 50,
       units = "cm",
       dpi=300)












for (i in pathways)
{
    genes <- A2.path[A2.path$Description == i, ]$core_enrichment
    print(genes)
#    heat <- t(scale(t(data_Promoters[genes,
#                                     annotations_A2$Sample_Name])))
    heat <- data_Promoters[genes,
                                     annotations_A2$Sample_Name]
    color <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)    
    pheatmap::pheatmap(heat,
                       treeheight_row      = 0,
                       annotation_col      = annotations_A2["Sample_Group"],
                       custering_distance_cols= "correlation",
                       custering_distance_rows= "correlation",
                       clustering_method  = "complete",
                       scale               = "row",
                       cutree_cols         = 2,
                       color               = color,
                       main                = i,
                       width               = 8,
                       height              = 6,
                       filename            = file.path(Origin, "Classifier",
                                                       paste0("A2_pathways_", i, ".png")))
}





pathways <- read.csv(file.path(Origin, "Interpretation/Kegg/Gse/Nbr_KeggGse_Promoters_according_comparisons.csv"))
A2Inactive.path <- pathways[pathways$comparison == "A2Inactive", ]
A2Inactive.path <- ParseEnrichResults(A2Inactive.path, column = "core_enrichment")
pathways <- c("Amyotrophic lateral sclerosis", "Hepatocellular carcinoma",
              "Huntington disease", "Prion disease", "Ribosome",
              "Shigellosis")

for (i in pathways)
{
    genes <- A2Inactive.path[A2Inactive.path$Description == i, ]$core_enrichment
    print(genes)
#    heat <- t(scale(t(data_Promoters[genes,
#                                     annotations_A2Inactive$Sample_Name])))
    heat <- data_Promoters[genes,
                                     annotations_A2Inactive$Sample_Name]
    color <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)    
    pheatmap::pheatmap(heat,
                       treeheight_row      = 0,
                       annotation_col      = annotations_A2Inactive["Sample_Group"],
                       custering_distance_cols= "correlation",
                       custering_distance_rows= "correlation",
                       clustering_method  = "complete",
                       scale               = "row",
                       cutree_cols         = 2,
                       color               = color,
                       main                = i,
                       width               = 8,
                       height              = 6,
                       filename            = file.path(Origin, "Classifier",
                                                       paste0("A2Inactive_pathways_", i, ".png")))
}










umap_func(heat,
          sampleTable = annotations_A2[c("Sample_Group", "Sample_Name")],
          highlightedVarName = "Sample_Group")




# create logistic regression model
logistic_model <- glm(var1 ~ var2,
                      data=data_Promoters[genes, annotations_A2$Sample_Name],
                      family=binomial)

#Data frame with hp in ascending order
Predicted_data <- data.frame(var2=seq(
                                      min(df$var2), max(df$var2),len=500))

# Fill predicted values using regression model
Predicted <- data$var1 = predict(
                   logistic_model, Predicted_data, type="response")

# Plot Predicted data and original data points
plot(var1 ~ var2, data=df)
lines(var1 ~ var2, Predicted_data, lwd=2, col="green")










library("readxl")
features <- read_excel(, "features_selected.xlsx"))


heat <- t(scale(t(M.Promoters.annotated[na.omit(features[['Logistic Regression']]),
                                        annotations_A2$Sample_Name])))
pheatmap::pheatmap(heat, treeheight_row = 0,
                   annotation_col = annotations_A2["Sample_Group"],
                   custering_distance = "euclidean",
                   clustering_method = "ward.D2",
                   cutree_cols = 2)


heat <- t(scale(t(M.Promoters.annotated[na.omit(features[['Random Forest']]),
                                        annotations_A2$Sample_Name])))
pheatmap::pheatmap(heat, treeheight_row = 0,
                   annotation_col = annotations_A2["Sample_Group"],
                   custering_distance = "euclidean",
                   clustering_method = "ward.D2",
                   cutree_cols = 2)

heat <- t(scale(t(M.Promoters.annotated[na.omit(features[['XGBoost']]),
                                        annotations_A2$Sample_Name])))
pheatmap::pheatmap(heat, treeheight_row = 0,
                   annotation_col = annotations_A2["Sample_Group"],
                   custering_distance = "euclidean",
                   clustering_method = "ward.D2",
                   cutree_cols = 2)



# all promoters DF in A2
heat <- t(scale(t(data_Promoters[pro.A2$Gene,
                                        annotations_A2$Sample_Name])))
pheatmap::pheatmap(heat, treeheight_row = 0,
                   annotation_col = annotations_A2["Sample_Group"],
                   custering_distance = "euclidean",
                   clustering_method = "ward.D2",
                   cutree_cols = 2)


na.omit(features[['Random Forest']])[ ! na.omit(features[['Random Forest']]) %in% pro.A2$Gene]

na.omit(features[['Logistic Regression']])[ ! na.omit(features[['Logistic Regression']]) %in% pro.A2$Gene]

na.omit(features[['XGBoost']])[ ! na.omit(features[['XGBoost']]) %in% pro.A2$Gene]




library("KEGGREST")
library("dplyr")
library("org.Hs.eg.db")

# Retrieve all kegg pathways, genes and pathway description
hsa_path_eg <- keggLink("pathway", "hsa") %>%
         tibble(pathway = ., eg = sub("hsa:", "", names(.)))

hsa_kegg_anno <- hsa_path_eg %>%
    mutate(
        SYMBOL = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID")#,
        #ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
    )

hsa_pathways <- keggList("pathway", "hsa") %>%
        tibble(pathway = names(.), description = .)
hsa_pathways$pathway <- paste0("path:", hsa_pathways$pathway)

FinalData <- left_join(hsa_kegg_anno, hsa_pathways)
FinalData <- FinalData %>% arrange(pathway, SYMBOL)  %>% dplyr::select(-eg)

ssData <- split(data.frame(FinalData)$SYMBOL, data.frame(FinalData)$description)



M.Promoters.annotated <- read.csv(file.path(DirMatrices, "M_Promoters_annotated.csv"), row.names = 1)
print("matrices of M values loaded")


library("GSVA")

ssGSEA <- gsva(as.matrix(data_Promoters[annotations_A2$Sample_Name]),
               ssData,
               method = "ssgsea",
               verbose = TRUE)


sel <- annotations_A2[(annotations_A2$Sample_Group == "ATM" & annotations_A2$Inactive == "inactive ATM") | annotations_A2$Sample_Group == "Non.ATM" , ]
    
data1 <- ssGSEA[, sel$Sample_Name]
design <- model.matrix(~ factor(sel$Sample_Group, levels=c("ATM", "Non.ATM")))
colnames(design) <- c("Intercept", "ATMVsNon.ATM")
fit <- limma::lmFit(data1, design)
fit <- limma::eBayes(fit)
sigPathways <- limma::topTable(fit, coef="ATMVsNon.ATM", number=Inf, p.value=0.05,
                               adjust="BH")
sigPathways <- sigPathways[sigPathways$adj.P.Val<0.0001,]
res <- limma::decideTests(fit, p.value=0.05)

heat <- t(scale(t(data1[rownames(sigPathways), ])))

pheatmap::pheatmap(heat, treeheight_row = 0,
                   annotation_col = annotations_A2["Sample_Group"],
                   custering_distance = "euclidean",
                   clustering_method = "ward.D2",
                   cutree_cols = 2)



data1 <- ssGSEA[,annotations_A2$Sample_Name]
design <- model.matrix(~ factor(annotations_A2$Sample_Group, levels=c("ATM", "Non.ATM")))
colnames(design) <- c("Intercept", "ATMVsNon.ATM")
fit <- limma::lmFit(data1, design)
fit <- limma::eBayes(fit)
sigPathways <- limma::topTable(fit, coef="ATMVsNon.ATM", number=Inf, p.value=0.05,
                               adjust="BH")
sigPathways <- sigPathways[sigPathways$adj.P.Val<0.0001,]
res <- limma::decideTests(fit, p.value=0.05)

heat <- t(scale(t(data1[rownames(sigPathways), ])))

pheatmap::pheatmap(heat, treeheight_row = 0,
                   annotation_col = annotations_A2[c("Sample_Group")],#, "Study", "ER_status",
                                                   #"PR_status", "subtype", "Inactive")],
                   custering_distance = "euclidean",
                   clustering_method = "ward.D2",
                   cutree_cols = 4)







pathways <- read.csv(file.path(Origin, "Interpretation/Kegg/Gse/Nbr_KeggGse_Promoters_according_comparisons.csv"))

P <- pathways[pathways$comparison == "A2Inactive", ]$Description

ssGSEA <- data.frame(ssGSEA)
rownames(ssGSEA) <- gsub(" - Homo sapiens \\(human\\)", "", rownames(ssGSEA))

H1 <- preComputePheatmap(data=ssGSEA[P,],
                         annotations = annotations_A2,
                         method_dist = "pearson",
                         method_clust = "ward.D2",
                         cutree_cols        = 2,
                         cutree_rows        = 2,
                         showRownames = T,
                         colAnno = c("Grade","subtype", "Supposed biallelic inactivation", "Sample_Group"),
                   )










ssGSVA <- gsva(as.matrix(M.Promoters.annotated),
               ssData, #[grepl("^[Homolo]", names(ssData))],
               method = "gsva",
               verbose = TRUE)


sel <- annotations_A2[(annotations_A2$Sample_Group == "ATM" & annotations_A2$Inactive == "inactive ATM") | annotations_A2$Sample_Group == "Non.ATM" , ]

rownames(ssGSVA) <- gsub(" - Homo sapiens \\(human\\)", "", rownames(ssGSVA))

data1 <- ssGSVA[, #A2.Pro.KeggGseNbr$Description,
                annotations_A2$Sample_Name]
design <- model.matrix(~ factor(annotations_A2$Sample_Group, levels=c("ATM", "Non.ATM")))
colnames(design) <- c("Intercept", "ATMVsNon.ATM")
fit <- limma::lmFit(data1, design)
fit <- limma::eBayes(fit)
sigPathways <- limma::topTable(fit, coef="ATMVsNon.ATM", number=Inf, p.value=0.05,
                               adjust="BH")
sigPathways <- sigPathways[sigPathways$adj.P.Val<0.01,]
res <- limma::decideTests(fit, p.value=0.05)

heat <- t(scale(t(data1[rownames(sigPathways), ])))

pheatmap::pheatmap(heat, treeheight_row = 0,
                   annotation_col = annotations_A2["Sample_Group"],
                   custering_distance = "euclidean",
#                   clustering_distance_cols = "correlation",
                   clustering_method = "ward.D2",
                   cutree_cols = 2)








#heat <- t(scale(t(ssGSVA)))


test <- pheatmap::pheatmap(heat[c("Homologous recombination - Homo sapiens (human)",
                                  "Non-homologous end-joining - Homo sapiens (human)",
                                  "Viral carcinogenesis - Homo sapiens (human)"), ],
                           clustering_distance_cols = "correlation",
                           annotation_col = annotations_A2["Sample_Group"])


ggsave(file.path(OutFiles1, "test_ssGSVA.png"), plot = test)



MethylResolver <- function(){
    
library(MethylResolver)

load(file.path(Origin, "matrices/MethylSet_Funnorm_Normalisation.RData"))

load(file.path(Origin, "matrices/Beta_Funnorm_Normalisation.RData"))
#all_Beta <- read.csv(file.path("matrices/Beta_Funnorm_Normalisation.RData"), row.names=1)


# Deconvolution with default signature and only calculating relative fractions:
MethylResolver(methylMix = all_beta[,c(1:19, 21:26, 28:251, 253:260, 450:520)], absolute = FALSE, betaPrime = FALSE)

#1:19, 21:26, 28:251, 253:260
#480:520
deconv <- read.table(file = file.path(getwd(), "MethylResolver.txt"), sep="\t", header=TRUE)
deconv$Sample <- rownames(deconv)
deconv <- gather(deconv, CellType, score, colnames(deconv)[1]:colnames(deconv)[length(colnames(deconv)) - 1], factor_key=TRUE)
deconv <- merge(deconv, annotations["Sample_Group"], by.x="Sample", by.y="row.names")


# Relative fraction and Tumour Purity-Adjusted Fraction
ggplot(deconv[! deconv$CellType %in% c("RMSE1", "R1", "RMSE2", "R2"),],
       aes(fill=CellType, y=score, x=Sample)) +
    geom_bar(position="fill", stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

}

#dotplot(A1.Pro.KeggGse, showCategory=45) + scale_y_discrete(labels=function(x) str_wrap(x, width=150))

# working to compare hyper and hypomethylated pathways
#dotplot(A1.Pro.KeggGse, showCategory=45, split = ".sign") + facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
# Nouveaux noms d'étiquettes des facettes pour la variable `dose`
#slot(A1.Pro.KeggGse, "result")$NES_sign <- ifelse(slot(A1.Pro.KeggGse, "result")$NES < 0, 'hypomethylated', "hypermethylated")
#dotplot(A1.Pro.KeggGse, showCategory=45, split = "NES_sign") + facet_grid(.~NES_sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=50))




# Comparison of A1 and A2 (based on A1)
test <- A1.Pro.KeggGse

IDs <- slot(test, "result")$Description[1:45]

slot(test, "result") <- rbind(slot(test, "result")[1:45,], subset(slot(A2.Pro.KeggGse, "result"), slot(A2.Pro.KeggGse, "result")$Description %in% IDs))

slot(test, "result")$length_core <- str_count(test$core_enrichment, pattern = "/") + 1

slot(test, "result")$GeneRatio <- slot(test, "result")$length_core / slot(test, "result")$setSize

Order <- slot(test, "result")[order(slot(test, "result")$GeneRatio), ]
Order <- Order[Order$comparison == "A1", ]$Description

dot <- dotplot(test, showCategory=45, split="comparison") + facet_grid(.~comparison) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=150), limits = Order)

#dot

ggsaveDotPlot(filename = paste0(OutFiles, "GSE_comparison_A1_A2.png"),
              plot = dot)


# Comparison of A1 and A2 (based on both)
test <- A1.Pro.KeggGse

IDs <- slot(test, "result")$Description[1:45]

slot(test, "result") <- rbind(slot(test, "result")[1:45,], subset(slot(A2.Pro.KeggGse, "result"), slot(A2.Pro.KeggGse, "result")$Description %in% IDs))

slot(test, "result")$length_core <- str_count(test$core_enrichment, pattern = "/") + 1

slot(test, "result")$GeneRatio <- slot(test, "result")$length_core / slot(test, "result")$setSize

Order <- slot(test, "result")[order(slot(test, "result")$GeneRatio), ]
Order <- Order[Order$comparison == "A1", ]$Description

dot <- dotplot(test, showCategory=45, split="comparison") + facet_grid(.~comparison) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=150), limits = Order)

#dot

ggsaveDotPlot(filename = paste0(OutFiles, "GSE_comparison_A1_A2.png"),
              plot = dot)



a <- dotplot(A1.Pro.KeggGse, showCategory=45, split = "NES_sign") + scale_y_discrete(labels=function(x) str_wrap(x, width=50)) + ggtitle("A1")

b <- dotplot(A1.Pro.KeggGse, showCategory=45, split = "NES_sign") + scale_y_discrete(labels=function(x) str_wrap(x, width=50)) + ggtitle("A2")

ggpubr::ggarrange(a, b, ncol=2)


a <- clusterProfiler::dotplot(A1.Pro.KeggGse, showCategory=40) + scale_y_discrete(labels=function(x) str_wrap(x, width=50)) + ggtitle("A1")

b <- enrichplot::dotplot(A2.Pro.KeggGse, showCategory=40) + scale_y_discrete(labels=function(x) str_wrap(x, width=50)) + ggtitle("A2")

ggpubr::ggarrange(a, b, ncol=2)



load(file.path(OutFiles, "AllKeggGseObjects.Rdata"))

getGeneRatio <- function(object)
    {
        ratio <- c()
        for (i in 1:length(object$"Description")){
            ratio <- c(ratio, (stringr::str_count(object$core_enrichment[i], "/") + 1) / object$setSize[i])
        }
        slot(object, "result")$GeneRatio <- ratio
        return(object)
    }


A1.Pro.KeggGse <- getGeneRatio(A1.Pro.KeggGse)

A2.Pro.KeggGse <- getGeneRatio(A2.Pro.KeggGse)

Res <- as.data.frame(A1.Pro.KeggGse)
idx <- order(Res$GeneRatio, decreasing=TRUE)
Res1 <- Res[idx,]

Res <- as.data.frame(A2.Pro.KeggGse)
idx <- order(Res$GeneRatio, decreasing=TRUE)
Res2 <- Res[idx,]


clusterProfiler::dotplot(A1.Pro.KeggGse, showCategory=45)

dotplot_custom(as.data.frame(A1.Pro.KeggGse), showCategory=50)

dotplot_custom(as.data.frame(A2.Pro.KeggGse), showCategory=50)


# Specific to A1 (ER,)
match(Res1$Description, Res2$Description)
Res1$Description[which(match(Res1$Description, Res2$Description) %in% NA)]


# Specific to A2 (ATM tumours)
match(Res2$Description, Res1$Description)
Res2$Description[which(match(Res2$Description, Res1$Description) %in% NA)]



dim(as.data.frame(A1.Pro.KeggGse))








#plot <- lolliplot_custom_Paper(data = rbind(as.data.frame(A1.Pro.KeggGse),
#                                            as.data.frame(A2.Pro.KeggGse)
#                                            ),
#                       main              = NULL,
#                       comparison        = TRUE,
#                       orderAccordingTo  = "GeneRatio",
#                       size              = "NES",
#                       showSign          = T,
#                       showCategory      = NULL,
#                       sign              = NULL)

#plot <- plot + scale_x_discrete(labels = function(x) str_wrap(x, width = 40))
#plot <- plot + scale_y_discrete(labels = function(y) str_wrap(y, width = 60))
#ggsave(file.path(OutFiles, "lolliplot_comparison.png"), plot, width = 13.2, height=15, units="in")






eeeeeeeee








# Create files for Cytoscape at the pathway level

file <- file.path(Origin, "Interpretation/Kegg/PathwayOver-representation/Kegg_Promoters_according_comparisons.csv")
write.gmt(file, type = "ORA", level = "pathway", sep="\t")

#file <- file.path(Origin, "Interpretation/Kegg/PathwayOver-representation/Kegg_Genes_according_comparisons.csv")
#write.gmt(file, type = "ORA", level = "pathway", sep=",")


file <- file.path(Origin, "/Interpretation/Kegg/Gse/KeggGse_Promoters_according_comparisons.csv")
write.gmt(file, type = "GSEA", level = "pathway", sep=",")

#file <- file.path(Origin, "/Interpretation/Kegg/Gse/KeggGse_Genes_according_comparisons.csv")
#write.gmt(file, type = "GSEA", level = "pathway", sep="\t")

file <- file.path(Origin, "Interpretation/Kegg/Gse/Nbr_KeggGse_Promoters_according_comparisons.csv")
write.gmt(file, type = "GSEA", level = "pathway", sep=",")


file <- file.path(Origin, "Interpretation/Kegg/Gse/KeggGseFilter_Promoters_according_comparisons.csv")
write.gmt(file, type = "GSEA", level = "pathway", sep=",")


#file <- "~/ATM_Analysis/svn/analyse/results/FinalMethylation/Interpretation/Kegg/Gse/KeggGse_Promoters_AOld.csv"
#write.gmt(file, type = "GSEA", level = "pathway")

#file <- "~/ATM_Analysis/svn/analyse/results/FinalMethylation/Interpretation/Kegg/Gse/KeggGse_Promoters_AOldER.csv"
#write.gmt(file, type = "GSEA", level = "pathway")



# Create files for Cytoscape at the gene level
file <- file.path(Origin, "Interpretation/Kegg/PathwayOver-representation/Kegg_Promoters_according_comparisons.csv")
write.gmt(file, type = "ORA", level = "gene", sep="\t")

#file <- file.path(Origin, "Interpretation/Kegg/PathwayOver-representation/Kegg_Genes_according_comparisons.csv")
#write.gmt(file, type = "ORA", level = "gene", sep=",")


file <- file.path(Origin, "Interpretation/Kegg/Gse/KeggGse_Promoters_according_comparisons.csv")
write.gmt(file, type = "GSEA", level = "gene", sep=",")

#file <- file.path(Origin, "Interpretation/Kegg/Gse/KeggGse_Genes_according_comparisons.csv")
#write.gmt(file, type = "GSEA", level = "gene", sep="\t")


file <- file.path(Origin, "Interpretation/Kegg/Gse/Nbr_KeggGse_Promoters_according_comparisons.csv")
write.gmt(file, type = "GSEA", level = "gene", sep=",")



stop("enrichemnt performed")






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
            print(EnrichmentObject$Description[pathN])
            # GSEA plot
            gseaN <- gseaplot2(EnrichmentObject,
                               geneSetID     = pathN,
                               title         = EnrichmentObject$Description[pathN],
                               pvalue_table  = TRUE)

            print("save of plot 1")
            print(EnrichmentObject$Description[pathN])
            ggsave(file.path(OutFiles,
                             paste0(name, "_GseaPlot_", EnrichmentObject$Description[pathN], '.png')),
                   gseaN,
                   height  = 7,
                   width   = 14)
            # Get the list of genes for the enriched terms
            List <- ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description == EnrichmentObject$Description[pathN],],
                                       column = "core_enrichment")
            print("list 1")
            print(List)
            List <- List$core_enrichment
            print("list 2")
            print(head(List))
            
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
            print(List)
            heatmapN <- heatmap_func(data          = dataFrame.annotated[List,],
                                     annotations   = annotations,
                                     fontsize_row  = 8,
                                     fontsize_col  = 10,
                                     color         = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),

                                     main          = EnrichmentObject$Description[pathN],
                                     filename      = file.path(OutFiles,
                                                               paste0(name, "_M_Heatmap_", EnrichmentObject$Description[pathN], '.png')))

#            print("save of plot 2")
#            return(heatmapN)
#            ggsave(file.path(OutFiles,
#                             paste0(name, "_M_Heatmap_", EnrichmentObject$Description[pathN], '.png')),
#                   plot = heatmapN)#,
                   #height  = 2100,
                   #width   = 4200,
                   #units = "px")

            heatmapN2 <- heatmap_func(data         = dataFrame.annotated2[List,],
                                     annotations   = annotations,
                                     fontsize_row  = 8,
                                     fontsize_col  = 10,
                                     color         = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),

                                     main          = EnrichmentObject$Description[pathN],
                                     filename      = file.path(OutFiles,
                             paste0(name, "_Beta_Heatmap_", EnrichmentObject$Description[pathN], '.png')))

           # print("save of plot 3")
           # ggsave(,
           #        heatmapN2)#,
                   #height  = 2100,
                   #width   = 4200,
                   #units = "px")
            
            # UMAP of corresponding genes
#            k <- 2
#            umap.plot <- umap_func(matrix  = dataFrame.annotated[List,],
#                                   sampleTable         = annotations,
#                                   highlightedVar      = "Sample_Group",
#                                   highlightedVarName  = "Sample_Group",
#                                   colors              = c("orange", "Red"),
#                                   sampleLabel         = "Sample_Name",
#                                   metric              = "euclidean",
#                                   title               = EnrichmentObject$Description[pathN],
#                                   fontsize            = 14,
#                                   assignGroup         = TRUE,
#                                   k = k
#                                   )

 #           print("save of plot 4")
 #           ggsave(filename = file.path(OutFiles,
 #                                       paste0(name, "_UMAP_", EnrichmentObject$Description[pathN], '.png')),
 #                  plot = umap.plot,
 #                  width = 12, height = 7)
        }
    }
}



# Annotations
#annotations <- read.csv(file.path("/data/users/nviart/ATM_Analysis/svn/analyse/results/methylation/matrices/", "Annotations.csv"), row.names = 1)

OutFiles = file.path(OutFiles1, "Kegg/Gse/")
load(file.path(OutFiles, "AllKeggGseObjects.Rdata"))

# Recover raw data to plot heatmaps
# Matrices of M values annotated for genes
M.Genes.annotated <- read.csv(file.path(DirMatrices, "M_Genes_annotated.csv"), row.names = 1)

# Matrices of M values annotated for promoters
M.Promoters.annotated <- read.csv(file.path(DirMatrices, "M_Promoters_annotated.csv"), row.names = 1)
print("matrices of M values loaded")

# Matrices of M values annotated for genes
Beta.Genes.annotated <- read.csv(file.path(DirMatrices, "Beta_Genes_annotated.csv"), row.names = 1)

# Matrices of M values annotated for promoters
Beta.Promoters.annotated <- read.csv(file.path(DirMatrices, "Beta_Promoters_annotated.csv"), row.names = 1)
print("matrices of Beta values loaded")


# Plot all individual enrichment terms (GSEA, heatmap and umap)
OutFiles = file.path(OutFiles1, "Kegg/Gse/IndividualPlots/")

if(!dir.exists(OutFiles))
{
    dir.create(OutFiles, recursive = TRUE)
}

print("beginning of first GSEA")
plotGSEA(c("A2.Pro.KeggGse", "A0A2.Pro.KeggGse", "A2Inactive.Pro.KeggGse"))

print("first GSEA done")
plotGSEA(c("A2.Genes.KeggGse", "A0A2.Genes.KeggGse", "A2Inactive.Genes.KeggGse"))
print("second GSEA done")












































# Plot of several enrichemnt terms for the same Analysis

####### For a particular analysis and several pathways

AggregatePathways <- function(EnrichmentObject,
                              ListOfPathways,
                              type     = "all",
                              freq     = NULL,
                              k        = 2,
                              samples  = NULL # By default all samples are used
)
{
    # get informations for the name of the file
    Feature <- strsplit(EnrichmentObject, "\\.")[[1]][2]
    Analysis <- strsplit(EnrichmentObject, "\\.")[[1]][1]
    name <- paste0(Analysis, ".", Feature)

    # Evaluate the object (from string to object)
    EnrichmentObject <- eval(parse(text = EnrichmentObject))
    # Verification of pathways found
    print(EnrichmentObject[EnrichmentObject$Description %in% ListOfPathways,]$Description)
    # Get the gene names for each pathways indicated as input
    List <- ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description %in% ListOfPathways,],
                               column = "core_enrichment")
    # Handle if you want to filter genes to have genes present in a certain percentage of pathways
    if (type == "intercept")
    {
        Freq <- table(List$core_enrichment)
        if (! is.null(freq))
        {
            Genes <- names(Freq[Freq >= length(ListOfPathways)*freq])
        } else {
            message("Be CAREFULL, freq parameter is not defined and all genes of all found pathways will be used")
            Genes <- List$core_enrichment
        }
        if (length(Genes) == 0)
        {
            stop("Your criteria is too much stringent, no genes are retained and no Heatmap will be created")
        }
    } else {
        Genes <- List$core_enrichment
    }
    # Suppress duplicated genes
    Genes <- unique(Genes)
    print(length(Genes))
    
    if (dim(List)[1] != 0){
        # Handle if the filename is too long
        fileName = paste0(unique(List$Description), collapse = "\t")
        fileName = stringr::str_replace_all(fileName, " - ", "AND")
        fileName = stringr::str_replace_all(fileName, "-", "_")
        fileName = stringr::str_replace_all(fileName, " ", "_")
        if (nchar(fileName) > 20)
        {
            fileName <- paste0(substr(fileName, 1, 60))
        }
        # Which Feature is used ==> choose the correct dataset
        if (Feature == "Pro")
        {
            dataset = M.Promoters.annotated
        }
        if (Feature == "Genes")
        {
            dataset = M.Genes.annotated
        }
        # handle if you want only a subset of samples
        if (! is.null(samples))
        {
            #print(samples)
            dataset = dataset[, samples]
            annotations = subset(annotations, annotations$Sample_Name %in% samples)
        }
        # Test with manual clustering
        #dataset2 <- dataset[Genes,]
        #dataset2 <- t(dataset2[complete.cases(dataset2),])
        #dataset2 <- scale(dataset2)
        #print(colnames(dataset2))
        #d <- dist(dataset2, method = "euclidean")
        #png(file.path(OutFiles,
        #              paste0(name,
        #                     "_Dendo_",
        #                     ifelse(!is.null(freq), freq, ""),
        #                     "_", 
        #                     fileName,
        #                     ".png")
        #              ),
        #    width = 700, height = 600
        #    )
        #hc <- hclust(d, method = "complete")
        #plot(hc, hang = -1)
        #dev.off()

        #print(colnames(dataset))
        #print(annotations$Sample_Name)
        # Create and save the heatmap
        main = paste0(unique(List$Description), collapse = " & ")
        heatmap <- heatmap_func(data              = dataset[Genes,],
                                annotations       = annotations,
                                subsetAnnotation  = c("Study", "ER_status", "Variant_type", "ATM_LOH"),
                                color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                                fontsize_row      = 8,
                                fontsize_col      = 10,
                                main              = main,
                                filename          = file.path(OutFiles,
                                    paste0(name,
                                           "_Heatmap_",
                                           ifelse(!is.null(freq), freq, ""),
                                           "_",
                                           ifelse(!is.null(samples),
                                                  stringr::str_split(
                                                      stringr::str_split(
                                                          deparse(substitute(annotations_A1$Sample_Name)), "_")[[1]][2], "\\$")[[1]][1]
                                                  
                                                  
                                                , ""),
                                           fileName,
                                           ".png")
                                                     )
                                )
        # Save the list of genes
        listGenes <- as.data.frame(sort(Genes[Genes %in% rownames(dataset)])) 
        write.table(listGenes,
                    file = file.path(OutFiles,
                                    paste0(name,
                                           "_Genes_",
                                           ifelse(!is.null(freq), freq, ""),
                                           "_",
                                           ifelse(!is.null(samples),
                                                  stringr::str_split(
                                                      stringr::str_split(
                                                          deparse(substitute(annotations_A1$Sample_Name)), "_")[[1]][2], "\\$")[[1]][1]
                                                  
                                                  
                                                , ""),
                                           fileName,
                                           ".csv")
                                                     ),
                    quote = FALSE,
                    sep = ",",
                    row.names = FALSE,
                    col.names = FALSE)
        #main = paste0(unique(List$Description), collapse = " & ")
        #heatmap <- heatmap_func(data          = dataset[Genes,],
        #                        annotations   = annotations,
        #                        fontsize_row  = 8,
        #                        fontsize_col  = 10,
        #                        main = main)

        #k <- 3
        #clusters <- cutree(heatmap$tree_col, k=3)

        #target <- ifelse(annotations$Sample_Group == "ATM", 1, 2)
        #names(target) <- annotations$Sample_Name

        #confusionM <- table(clusters, target)
        #clusterNoControl <- names(confusionM[,2][confusionM[,2] == 0])
        #clusterToAdd <- names(confusionM[,1][confusionM[,1] == max(confusionM[,1])])

        #if (k == 3)
        #{
        #    clusters[clusters == clusterNoControl] <- clusterToAdd
        #}

        #cm <- table(clusters, target)
        #print("Rand Index")
        #print(sum(diag(cm)) / sum(cm))
        #print("Normalised Mutual Information")
        #print(aricode::NMI(clusters, target, variant = "max"))
        
        #ggsave(filename = file.path(OutFiles,
        #                            paste0(name,
        #                                   "_Heatmap2_",
        #                                   ifelse(!is.null(freq), freq, ""),
        #                                   "_", 
        #                                   fileName,
        #                                   ".png")
        #                            ),
        #       plot = heatmap
        #       )

        print("heatmap done and saved")
        # Create and save the umap
        k <- k
        #print("test1")
        #print(paste(unique(List$Description), collapse = " & "))

        #umap.plot <- umap_func(matrix  = dataset[Genes,],
        #                       sampleTable         = annotations,
        #                       highlightedVar      = "Sample_Group",
        #                       highlightedVarName  = "Sample_Group",
        #                       colors              = c("orange", "Red"),
        #                       sampleLabel         = "Sample_Name",
        #                       metric              = "euclidean",
        #                       title               = paste0(unique(List$Description), collapse = " & "),
        #                       fontsize            = 14,
        #                       assignGroup         = TRUE,
        #                       k                   = k
        #                       )
        #print("test2")
        #ggsave(filename = file.path(OutFiles,
        #                            paste0(name,
        #                                   "_Umap_",
        #                                   ifelse(!is.null(freq), freq, ""),
        #                                   "_",
        #                                   fileName,
        #                                   ".png")
        #                            ),
        #       plot = umap.plot)
    }
}  


OutFiles = file.path(OutFiles1, "Kegg/Gse/CombinedPlots_withoutCytoscape/")

if(!dir.exists(OutFiles))
{
    dir.create(OutFiles, recursive = TRUE)
}


########################
###### A4 analysis #####
########################

## Cancer diseases
#"Neutrophil extracellular trap formation"
#"Pathways in cancer"
#"Systemic lupus erythematosus"
#"Transcriptional misregulation in cancer"
#"Viral carcinogenesis"

ListOfPathways <- c("Neutrophil extracellular trap formation", "Transcriptional misregulation in cancer")
AggregatePathways("A4.Pro.KeggGse", ListOfPathways)

ListOfPathways <- c("Neutrophil extracellular trap formation", "Systemic lupus erythematosus","Transcriptional misregulation in cancer")
AggregatePathways("A4.Pro.KeggGse", ListOfPathways)

ListOfPathways <- c("Neutrophil extracellular trap formation", "Pathways in cancer", "Transcriptional misregulation in cancer", "Viral carcinogenesis")
AggregatePathways("A4.Pro.KeggGse", ListOfPathways)


# Diseases of aggregation
#"Amyotrophic lateral sclerosis"
#"Huntington disease"
#"Parkinson disease"
#"Prion disease"

ListOfPathways <- c("Huntington disease", "Prion disease")
AggregatePathways("A4.Pro.KeggGse", ListOfPathways)

ListOfPathways <- c("Prion disease", "Parkinson disease", "Huntington disease")
AggregatePathways("A4.Pro.KeggGse", ListOfPathways)

ListOfPathways <- c("Prion disease", "Parkinson disease", "Huntington disease", "Amyotrophic lateral sclerosis")
AggregatePathways("A4.Pro.KeggGse", ListOfPathways)


########################
###### A3 analysis #####
########################

ListOfPathways <- c("Neutrophil extracellular trap formation", "Systemic lupus erythematosus", "Transcriptional misregulation in cancer", "Viral carcinogenesis")
AggregatePathways("A3.Pro.KeggGse", ListOfPathways)

#A3.Pro.KeggGse$Description


########################
###### A2 analysis #####
########################

ListOfPathways <- c("Fanconi anemia pathway")
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)

ListOfPathways <- c("Fanconi anemia pathway", "Transcriptional misregulation in cancer")
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)

ListOfPathways <- c("Fanconi anemia pathway", "Neutrophil extracellular trap formation")
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)

ListOfPathways <- c("Fanconi anemia pathway", "Neutrophil extracellular trap formation", "Transcriptional misregulation in cancer")
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)

ListOfPathways <- c("Fanconi anemia pathway", "Basal transcription factors", "Transcriptional misregulation in cancer")
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)

ListOfPathways <- c("Transcriptional misregulation in cancer", "Basal transcription factors", "Fanconi anemia pathway", "Prion disease", "Parkinson disease", "Huntington disease")       
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)


########################
###### A1 analysis #####
########################







#test.Genes <- getAllGenesForOneEnrichedPathway(c("A1.Genes.KeggGse", "A2.Genes.KeggGse", "A3.Genes.KeggGse", "A4.Genes.KeggGse"))


#List <- ParseEnrichResults(enrichResultObject = A1.Pro.KeggGse[A1.Pro.KeggGse$Description %in% c("p53 signaling pathway", "PI3K-Akt signaling pathway", "Transcriptional misregulation in cancer", "Fanconi anemia pathway")],
#                           column = "core_enrichment")

#List <- ParseEnrichResults(enrichResultObject = A1.Pro.KeggGse[A1.Pro.KeggGse$Description == "Transcriptional misregulation in cancer"],
#                           column = "core_enrichment")

#List <- ParseEnrichResults(enrichResultObject = A4.Pro.KeggGse[A4.Pro.KeggGse$Description == "Transcriptional misregulation in cancer",],
#                           column = "core_enrichment")

#List <- c(List ,  ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description %in% c("Basal transcription factors", "Fanconi anemia pathway"),],
#                           column = "core_enrichment")
#)


#List <- c(ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description %in% c("Basal transcription factors", "Fanconi anemia pathway"),],
#                           column = "core_enrichment")
#)

#List <- c(ParseEnrichResults(enrichResultObject = A4.Pro.KeggGse[A4.Pro.KeggGse$Description %in% c("Transcriptional misregulation in cancer", "Basal transcription factors", "Fanconi anemia pathway"),],
#                           column = "core_enrichment")
#)




####### According to similarity terms
OutFiles = file.path(OutFiles1, "Kegg/Gse/CombinedPlotsFinal/")

if(!dir.exists(OutFiles))
{
    dir.create(OutFiles, recursive = TRUE)
}


load(file.path(OutFiles1, "Kegg/Gse/AllKeggGseObjects.Rdata"))

firstup <- function(x) {
    x <- tolower(x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    return(x)
  }

############################# Promoters
#Load masterPathways
Masters <- read.csv(file.path(OutFiles1, "Kegg/Gse/A1_MasterPathways.csv"))

##### A1
# Cluster 1 with only A1 samples
#ListOfPathways <- unlist(lapply(Masters$MasterPathway1, firstup))
ListOfPathways <- c("Lipid and atherosclerosis", "Melanoma", "Hepatitis B", "Pathways in cancer", "Kaposi sarcoma-associated herpesvirus infection", "Gastric cancer", "Human papillomavirus infection", "Basal cell carcinoma", "Prolactin signaling pathway", "Human cytomegalovirus infection", "Small cell lung cancer", "Hepatocellular carcinoma", "Cellular senescence", "Non-small cell lung cancer", "Chronic myeloid leukemia", "Proteoglycans in cancer", "Insulin signaling pathway", "EGFR tyrosine kinase inhibitor resistance", "Glioma", "Breast cancer", "Endocrine resistance", "Epstein-Barr virus infection", "Wnt signaling pathway", "Signaling pathways regulating pluripotency of stem cells", "Melanogenesis", "ErbB signaling pathway")

AggregatePathways("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name) 
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3,
#                  samples = annotations_A1$Sample_Name)

# Cluster 1
AggregatePathways("A1.Pro.KeggGse", ListOfPathways) 
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3)

# Cluster2 with only A1 samples
#ListOfPathways <- unlist(lapply(Masters$MasterPathway2, firstup))
ListOfPathways <- c("Huntington disease", "Pathways of neurodegeneration - multiple diseases", "Amyotrophic lateral sclerosis", "Alzheimer disease", "Spinocerebellar ataxia", "Parkinson disease", "Oxidative phosphorylation", "Thermogenesis", "Non-alcoholic fatty liver disease", "Chemical carcinogenesis - reactive oxygen species")
AggregatePathways("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3,
#                  samples = annotations_A1$Sample_Name)

# Cluster2
AggregatePathways("A1.Pro.KeggGse", ListOfPathways)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3)

# Cluster 3 with only A1 samples
#ListOfPathways <- unlist(lapply(Masters$MasterPathway3, firstup))
ListOfPathways <- c("Retinol metabolism", "Caffeine metabolism", "Steroid hormone biosynthesis", "Chemical carcinogenesis - DNA adducts", "Drug metabolism - cytochrome p450")
AggregatePathways("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3,
#                  samples = annotations_A1$Sample_Name)

# Cluster 3
AggregatePathways("A1.Pro.KeggGse", ListOfPathways)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3)




# Cluster 4 with only A1 samples
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3,
#                  samples = annotations_A1$Sample_Name)

# Cluster 4
AggregatePathways("A1.Pro.KeggGse", ListOfPathways)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3)



# Cluster 5 with only A1 samples
ListOfPathways <- c("Nucleotide metabolism", "Purine metabolism", "Pyrimidine metabolism")
AggregatePathways("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2,
#                  samples = annotations_A1$Sample_Name)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3,
#                  samples = annotations_A1$Sample_Name)

# Cluster 5
AggregatePathways("A1.Pro.KeggGse", ListOfPathways)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2)
#AggregatePathways("A1.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3)


# Cluster 6 with only A1 samples
ListOfPathways <- c("Cell cycle", "Oocyte meiosis")
AggregatePathways("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name)
# Cluster 6
AggregatePathways("A1.Pro.KeggGse", ListOfPathways)








####### A2
# Cluster 1 with only A2 samples
ListOfPathways <- c("Lipid and atherosclerosis", "Melanoma", "Hepatitis B", "Pathways in cancer", "Kaposi sarcoma-associated herpesvirus infection", "Gastric cancer", "Human papillomavirus infection", "Basal cell carcinoma", "Human cytomegalovirus infection", "Small cell lung cancer", "Hepatocellular carcinoma", "Cellular senescence", "Non-small cell lung cancer", "Chronic myeloid leukemia", "Proteoglycans in cancer", "Insulin signaling pathway", "EGFR tyrosine kinase inhibitor resistance", "Glioma", "Breast cancer", "Endocrine resistance", "Epstein-Barr virus infection", "Wnt signaling pathway", "Signaling pathways regulating pluripotency of stem cells", "Melanogenesis", "Hippo signaling pathway", "mTOR signaling pathway")

AggregatePathways("A2.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A2$Sample_Name)

# Cluster 1 
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)


# Cluster2 with only A1 samples
ListOfPathways <- c("Huntington disease", "Pathways of neurodegeneration - multiple diseases", "Amyotrophic lateral sclerosis", "Alzheimer disease", "Spinocerebellar ataxia", "Parkinson disease", "Oxidative phosphorylation", "Thermogenesis", "Non-alcoholic fatty liver disease", "Chemical carcinogenesis - reactive oxygen species")
AggregatePathways("A2.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A2$Sample_Name)

# Cluster2
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)


# Cluster 3 with only A2 samples
ListOfPathways <- c("Retinol metabolism", "Metabolism of xenobiotics by cytochrome P450", "Caffeine metabolism", "Steroid hormone biosynthesis", "Chemical carcinogenesis - DNA adducts", "Drug metabolism - cytochrome p450")
AggregatePathways("A2.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A2$Sample_Name)

# Cluster 3
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)


# Cluster4 with only A2 samples
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("A2.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A2$Sample_Name)

# Cluster4
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)



# Cluster 5 with only A2 samples
ListOfPathways <- c("Nucleotide metabolism", "Pyrimidine metabolism")
AggregatePathways("A2.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A2$Sample_Name)

# Cluster 5
AggregatePathways("A2.Pro.KeggGse", ListOfPathways)












####### AOld
# Cluster 1 with only AOld samples
ListOfPathways <- c("Parkinson disease", "Diabetic cardiomyopathy", "Thermogenesis", "Non-alcoholic fatty liver disease", "Chemical carcinogenesis - reactive oxygen species", "Amyotrophic lateral sclerosis", "Huntington disease", "Prion disease", "Pathways of neurodegeneration - multiple diseases", "Alzheimer disease", "Oxidative phosphorylation")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways,
                  samples = annotations_AOld$Sample_Name)

# Cluster 1
ListOfPathways <- c("Parkinson disease", "Diabetic cardiomyopathy", "Thermogenesis", "Non-alcoholic fatty liver disease", "Chemical carcinogenesis - reactive oxygen species", "Amyotrophic lateral sclerosis", "Huntington disease", "Prion disease", "Pathways of neurodegeneration - multiple diseases", "Alzheimer disease", "Oxidative phosphorylation")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways)


# Cluster 2 with only AOld samples
ListOfPathways <- c("Calcium signaling pathway", "Circadian entrainment", "Morphine addiction", "Cholinergic synapse", "Glutamatergic synapse")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways,
                  samples = annotations_AOld$Sample_Name)

# Cluster 2
ListOfPathways <- c("Calcium signaling pathway", "Circadian entrainment", "Morphine addiction", "Cholinergic synapse", "Glutamatergic synapse")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways)


# Cluster 3 with only AOld samples
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways,
                  samples = annotations_AOld$Sample_Name)

# Cluster 3
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways)


# Cluster 4 with only AOld samples
ListOfPathways <- c("2-Oxocarboxylic acid metabolism", "Biosynthesis of amino acids", "Citrate cycle (TCA cycle)", "Carbon metabolism")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways,
                  samples = annotations_AOld$Sample_Name)

# Cluster 4
ListOfPathways <- c("2-Oxocarboxylic acid metabolism", "Biosynthesis of amino acids", "Citrate cycle (TCA cycle)", "Carbon metabolism")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways)


# Cluster 5 with only AOld samples
ListOfPathways <- c("Arrhythmogenic right ventricular cardiomyopathy", "Dilated cardiomyopathy", "Adrenergic signaling in cardiomyocytes")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways,
                  samples = annotations_AOld$Sample_Name)

# Cluster 4A
ListOfPathways <- c("Arrhythmogenic right ventricular cardiomyopathy", "Dilated cardiomyopathy", "Adrenergic signaling in cardiomyocytes")
AggregatePathways("AOld.Pro.KeggGse", ListOfPathways)






####### AOld + ER
# Cluster 1 with only AOldER samples
ListOfPathways <- c("Prion disease", "Amyotrophic lateral sclerosis", "Oxidative phosphorylation", "Huntington disease", "Alzheimer disease", "Pathways of neurodegeneration - multiple diseases", "Parkinson disease")
AggregatePathways("AOldER.Pro.KeggGse", ListOfPathways,
                  samples = annotations_AOldER$Sample_Name)

# Cluster 1
ListOfPathways <- c("Prion disease", "Amyotrophic lateral sclerosis", "Oxidative phosphorylation", "Huntington disease", "Alzheimer disease", "Pathways of neurodegeneration - multiple diseases", "Parkinson disease")
AggregatePathways("AOldER.Pro.KeggGse", ListOfPathways)


# Cluster 2 with only AOld samples
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("AOldER.Pro.KeggGse", ListOfPathways,
                  samples = annotations_AOldER$Sample_Name)

# Cluster 2 with only AOld samples
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("AOldER.Pro.KeggGse", ListOfPathways)













####### A3
# Cluster 1
ListOfPathways <- c("Parkinson disease",  "Huntington disease", "Oxidative phosphorylation", "Amyotrophic lateral sclerosis", "Pathways of neurodegeneration - multiple diseases")
AggregatePathways("A3.Pro.KeggGse", ListOfPathways)
AggregatePathways("A3.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A3.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A3.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A3.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3)

# Cluster2
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("A3.Pro.KeggGse", ListOfPathways)
AggregatePathways("A3.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A3.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A3.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A3.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3)



####### A4
# Cluster 1
ListOfPathways <- c("Parkinson disease", "Amyotrophic lateral sclerosis", "Pathways of neurodegeneration - multiple diseases", "Alzheimer disease", "Huntington disease", "Oxidative phosphorylation", "Prion disease")
AggregatePathways("A4.Pro.KeggGse", ListOfPathways)
AggregatePathways("A4.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A4.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A4.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A4.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3)

# Cluster2
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation", "Necroptosis")
AggregatePathways("A4.Pro.KeggGse", ListOfPathways)
AggregatePathways("A4.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A4.Pro.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A4.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A4.Pro.KeggGse", ListOfPathways, "intercept", freq = 1/3)
















############################# Genes
##### A1
# Cluster 1
ListOfPathways <- c("Thermogenesis", "Prion disease", "Non-alcoholic fatty liver disease", "Pathways of neurodegeneration - multiple diseases",  "Amyotrophic lateral sclerosis", "Parkinson disease", "Alzheimer disease", "Huntington disease", "Chemical carcinogenesis - reactive oxygen species", "Oxidative phosphorylation", )
AggregatePathways("A1.Genes.KeggGse", ListOfPathways)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/3)

# Cluster2
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("A1.Genes.KeggGse", ListOfPathways)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/3)

# Cluster 3
ListOfPathways <- c("Ribosome", "Coronavirus disease - COVID-19")
AggregatePathways("A1.Genes.KeggGse", ListOfPathways)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/3)
#AggregatePathways("A1.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/10)


####### A2
# Cluster 1
ListOfPathways <- c("Parkinson disease", "Thermogenesis", "Huntington disease", "Prion disease", "Chemical carcinogenesis - reactive oxygen species", "Alzheimer disease", "Amyotrophic lateral sclerosis", "Non-alcoholic fatty liver disease","Oxidative phosphorylation", "Pathways of neurodegeneration - multiple diseases")
AggregatePathways("A2.Genes.KeggGse", ListOfPathways)
AggregatePathways("A2.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A2.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A2.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A2.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/3)

# Cluster2
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("A2.Genes.KeggGse", ListOfPathways)
AggregatePathways("A2.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A2.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A2.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A2.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/3)


####### A3
# Cluster 1
ListOfPathways <- c("Thermogenesis", "Oxidative phosphorylation", "Huntington disease", "Non-alcoholic fatty liver disease")
AggregatePathways("A3.Genes.KeggGse", ListOfPathways)
AggregatePathways("A3.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A3.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A3.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/3)

# Cluster2
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathways("A3.Genes.KeggGse", ListOfPathways)
AggregatePathways("A3.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A3.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A3.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A3.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/3)



####### A4
# Cluster 1
ListOfPathways <- c("Neutrophil extracellular trap formation", "Systemic lupus erythematosus", "Alcoholism", "Viral carcinogenesis")
AggregatePathways("A4.Genes.KeggGse", ListOfPathways)
AggregatePathways("A4.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.9)
AggregatePathways("A4.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.7)
AggregatePathways("A4.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/2)
AggregatePathways("A4.Genes.KeggGse", ListOfPathways, "intercept", freq = 1/3)








ListOfPathways <- c("Neutrophil extracellular trap formation", "Systemic lupus erythematosus", "Alcoholism", "Viral carcinogenesis")
AggregatePathways2("A3.Genes.KeggGse", ListOfPathways, "intercept", freq = 0.9)



EnrichmentObject <- "A4.Pro.KeggGse"
EnrichmentObject <- print(eval(parse(text = EnrichmentObject)))
    # Get the gene names for each pathways indicated as input
List <- ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description %in% ListOfPathways,],
                           column = "core_enrichment")
dataset <- M.Genes.annotated[List$core_enrichment,]
dataset <- dataset[complete.cases(dataset),]

dataset <- scale(dataset)
dataset <- t(dataset)

print(dataset)
d <- dist(dataset, method = "euclidean")
        #print(M.Promoters.annotated[List$core_enrichment,])
hc <- hclust(d, method = "complete")
#plot(hc)









AggregatePathways <- function(EnrichmentObject,
                              ListOfPathways,
                              samples  = NULL, # By default, all samples are used
                              type     = "all",
                              freq     = NULL,
                              k        = 2)
{
    # get informations for the name of the file
    Feature <- strsplit(EnrichmentObject, "\\.")[[1]][2]
    Analysis <- strsplit(EnrichmentObject, "\\.")[[1]][1]
    name <- paste0(Analysis, ".", Feature)

    # Evaluate the object (from string to object)
    EnrichmentObject <- eval(parse(text = EnrichmentObject))
    # Verification of pathways found
    print(EnrichmentObject[EnrichmentObject$Description %in% ListOfPathways,]$Description)
    # Get the gene names for each pathways indicated as input
    List <- ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description %in% ListOfPathways,],
                               column = "core_enrichment")
    # Handle if you want to filter genes to have genes present in a certain percentage of pathways
    if (type == "intercept")
        {
            Freq <- table(List$core_enrichment)
            if (! is.null(freq))
                {
                    Genes <- names(Freq[Freq >= length(ListOfPathways)*freq])
                } else {
                    message("Be CAREFULL, freq parameter is not defined and all genes of all found pathways will be used")
                    Genes <- List$core_enrichment
                }
            if (length(Genes) == 0)
                {
                    stop("Your criteria is too much stringent, no genes are retained and no Heatmap will be created")
                }
        } else {
            Genes <- List$core_enrichment
        }
    # Suppress duplicated genes
    Genes <- unique(Genes)
    if (dim(List)[1] != 0){
        # Handle if the filename is too long
        fileName = paste0(unique(List$Description), collapse = "\t")
        fileName = stringr::str_replace_all(fileName, " - ", "AND")
        fileName = stringr::str_replace_all(fileName, "-", "_")
        fileName = stringr::str_replace_all(fileName, " ", "_")
        if (nchar(fileName) > 20)
        {
            fileName <- paste0(substr(fileName, 1, 60))
        }
        # Which Feature is used ==> choose the correct dataset
        if (Feature == "Pro")
        {
            dataset = M.Promoters.annotated
        }
        if (Feature == "Genes")
        {
            dataset = M.Genes.annotated
        }
        # handle if you want only a subset of samples
        if (! is.null(samples))
            {
                dataset = dataset[, samples]
                annotations = subset(annotations, annotations$Sample_Name %in% samples)
            }
        # Test with manual clustering
        dataset2 <- dataset[Genes,]
        dataset2 <- t(dataset2[complete.cases(dataset2),])
        dataset2 <- scale(dataset2)
        print(colnames(dataset2))
        d <- dist(dataset2, method = "euclidean")
        png(file.path(OutFiles,
                      paste0(name,
                             "_Dendo_",
                             ifelse(!is.null(freq), freq, ""),
                             "_", 
                             fileName,
                             ".png")
                      ),
            width = 700, height = 600
            )
        hc <- hclust(d, method = "complete")
        plot(hc, hang = -1)
        dev.off()

        # Create and save the heatmap
        main = paste0(unique(List$Description), collapse = " & ")
        heatmap <- heatmap_func(data          = dataset[Genes,],
                                annotations   = annotations,
                                fontsize_row  = 8,
                                fontsize_col  = 10,
                                main = main,
                                filename = file.path(OutFiles,
                                                     paste0(name,
                                                            "_Heatmap_",
                                                            ifelse(!is.null(freq), freq, ""),
                                                            "_", 
                                                            fileName,
                                                            ".png")
                                                         )
                                    )






        main = paste0(unique(List$Description), collapse = " & ")
        heatmap <- heatmap_func(data          = dataset[Genes,],
                                annotations   = annotations,
                                fontsize_row  = 8,
                                fontsize_col  = 10,
                                main = main)

        k <- 2
        clusters <- cutree(heatmap$tree_col, k=3)
        # define target clusters
        target <- ifelse(annotations$Sample_Group == "ATM", 2, 1)
        names(target) <- annotations$Sample_Name
        
        confusionM <- table(clusters, target)
        clusterNoControl <- names(confusionM[,2][confusionM[,2] == 0])
        clusterToAdd <- names(confusionM[,1][confusionM[,1] == max(confusionM[,1])])
        makeCluster <- function(hc, k)
            {
                clusters <- cutree(hc, k=k)
                confusionM <- table(clusters, target)
                return(confusionM)
            }
        
        confusionM <- makeCluster(hc,
                                  k = k)
        maxInit <- MAX <- max(confusionM)
        while (maxInit == MAX)
            {
                if(k<9){
                    confusionM <- makeCluster(hc,
                                              k = k + 1)
                    k <- k + 1
                    #print(confusionM)
                    MAX <- max(confusionM)
                    print(MAX)
                }
            }
        clusters <- cutree(hc, k=k)    
        confusionM <- table(clusters, target)
        # select clusters that have more than half curie samples and others
        rows <- rownames(confusionM)[confusionM[, 2] > (confusionM[, 2] + confusionM[, 1])/2]
        notRows <- rownames(confusionM)[! rownames(confusionM) %in% rows]
        # Replace data in clusters matrix
        clusterInRows <- clusters[clusters %in% rows]
        clusterInRows <- replace(clusterInRows, clusterInRows %in% rows, 2)
        clusterNotInRows <- clusters[clusters %in% notRows]
        clusterNotInRows <- replace(clusterNotInRows, clusterNotInRows %in% notRows, 1)
        clusters <- c(clusterInRows, clusterNotInRows)
        # New confusion matrix
        confusionM <- table(clusters, target)
        # Metrics
        RI <- sum(diag(confusionM)) / sum(confusionM)
        NMIscore <- aricode::NMI(clusters, target, variant = "sum")
        RImclust <- mclust::adjustedRandIndex(clusters, target)

        print("Rand Index")
        print(RI)
        print("Normalised Mutual Information")
        print(NMIscore)
        print("Rand Index mclust")
        print(RImclust)
        
        # color samples taken into account
        clustersCol <- clusters2[annotations$Sample_Name[heatmap$tree_col$order]]
        heatmap$gtable$grobs[[4]]$gp <- grid::gpar(col = clustersCol)

        #dend <- as.dendrogram(hc)
        #dendextend::labels_colors(dend) <- ""
        #clusters <- clusters[names(dendextend::labels_colors(dend))]
        #dendextend::labels_colors(dend) <- clusters

        ggsave(filename = file.path(OutFiles,
                                    paste0(name,
                                           "_Heatmap_",
                                           ifelse(!is.null(freq), freq, ""),
                                           "_", 
                                           fileName,
                                           ".png")
                                    ),
               plot = heatmap
               )

        
        print("heatmap done and saved")
        # Create and save the umap
        k <- k
        
        umap.plot <- umap_func(matrix  = dataset[Genes,],
                               sampleTable         = annotations,
                               highlightedVar      = "Sample_Group",
                               highlightedVarName  = "Sample_Group",
                               colors              = c("orange", "Red"),
                               sampleLabel         = "Sample_Name",
                               metric              = "euclidean",
                               title               = paste(unique(List$Description), collapse = " & "),
                               fontsize            = 14,
                               assignGroup         = TRUE,
                               k                   = k
                               )
        ggsave(filename = file.path(OutFiles,
                                    paste0(name,
                                           "_Umap_",
                                           ifelse(!is.null(freq), freq, ""),
                                           "_",
                                           fileName,
                                           ".png")
                                    ),
               umap.plot,
               width = 12, height = 7)
    }
}  




EnrichmentObject <- "A2.Pro.KeggGse"
EnrichmentObject <- print(eval(parse(text = EnrichmentObject)))
    # Get the gene names for each pathways indicated as input
List <- ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description %in% ListOfPathways,],
                           column = "core_enrichment")
dataset <- M.Genes.annotated[List$core_enrichment,]
dataset <- dataset[complete.cases(dataset),]

heatmap <- heatmap_func(data          = dataset,
                                annotations   = annotations,
                                fontsize_row  = 8,
                                fontsize_col  = 10,
                                main = ""
                                    )

        k <- 2
        clusters <- cutree(heatmap$tree_col, k=k)
        # define target clusters
        target <- ifelse(annotations$Sample_Group == "ATM", 2, 1)
        names(target) <- annotations$Sample_Name
        
        confusionM <- table(clusters, target)
        #clusterNoControl <- names(confusionM[,2][confusionM[,2] == 0])
        #clusterToAdd <- names(confusionM[,1][confusionM[,1] == max(confusionM[,1])])

makeCluster <- function(hc, k)
            {
                clusters <- cutree(hc, k=k)
                confusionM <- table(clusters, target)
                return(confusionM)
            }
        
        confusionM <- makeCluster(hc = heatmap$tree_col,
                                  k = k)
        maxInit <- MAX <- max(confusionM)
        while (maxInit == MAX)
            {
                if(k<9){
                    confusionM <- makeCluster(hc = heatmap$tree_col,
                                              k = k + 1)
                    k <- k + 1
                    print(confusionM)
                    MAX <- max(confusionM)
                    print(MAX)
                }
            }
        clusters <- cutree(tree = heatmap$tree_col,
                           k=k)    
        confusionM <- table(clusters, target)
        # select clusters that have more than half curie samples and others
        rows <- rownames(confusionM)[confusionM[, 2] > (confusionM[, 2] + confusionM[, 1])/2]
        notRows <- rownames(confusionM)[! rownames(confusionM) %in% rows]
        # Replace data in clusters matrix
        clusterInRows <- clusters[clusters %in% rows]
        clusterInRows <- replace(clusterInRows, clusterInRows %in% rows, 2)
        clusterNotInRows <- clusters[clusters %in% notRows]
        clusterNotInRows <- replace(clusterNotInRows, clusterNotInRows %in% notRows, 1)
        clusters2 <- c(clusterInRows, clusterNotInRows)
        # New confusion matrix
        confusionM <- table(clusters2, target)
        # Metrics
        RI <- sum(diag(confusionM)) / sum(confusionM)
        NMIscore <- aricode::NMI(clusters2, target, variant = "sum")
        RImclust <- mclust::adjustedRandIndex(clusters2, target)

        print("Rand Index")
        print(RI)
        print("Normalised Mutual Information")
        print(NMIscore)
        print("Rand Index mclust")
        print(RImclust)


# Modify colors of samples to indicate the ones used for calculation
clustersCol <- clusters2[annotations$Sample_Name[heatmap$tree_col$order]]
heatmap$gtable$grobs[[4]]$gp <- grid::gpar(col = clustersCol)
# Modify title
my_title <- paste0(unique(List$Description), collapse = " & ")
heatmap$gtable$grobs[[1]] <- grid::textGrob(my_title, hjust = 0.5)


my_text <- grid::textGrob(label = paste0("Rand Index: ", RI, " ; Normalised Mutual Information: ", NMIscore))



grid::grid.draw(gtable::gtable_add_grob(heatmap$gtable$grobs[[8]], my_text, 3, 2))




my_title <- textGrob(my_title)
pheatmap <- pheatmap(data_subset, silent = TRUE)
grid.arrange(grobs = list(my_title, one[[4]]), heights = c(0.1, 1))







heatmap


        # color samples taken into account
        dend <- as.dendrogram(heatmap$tree_col)
        dendextend::labels_colors(dend) <- ""
        clusters <- clusters[names(dendextend::labels_colors(dend))]
        dendextend::labels_colors(dend) <- clusters


plot(as.dendrogram(heatmap$tree_col))
rect.hclust(heatmap$tree_col, k=k
)



cols=dataset[order(match(rownames(dataset), pheatmap$gtable$grobs[[4]]$label)), ]$colors


clusters2 <- clusters[pheatmap$gtable$grobs[[4]]$label]

heatmap$gtable$grobs[[4]]$gp <- grid::gpar(col = clustersCol, fontsize = 20)


eeee




















d <- read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation/Interpretation/Kegg/Gse/KeggGse_Promoters_according_comparisons_Genes_A2.gmt", sep="\t", header=F)

unique(d["V1"])

genes <- subset(d, V1 %in% c("Chemical carcinogenesis - reactive oxygen species", "Prion disease", "Non-alcoholic fatty liver disease", "Amyotrophic lateral sclerosis", "Alzheimer disease", "Parkinson disease", "Pathways of neurodegeneration - multiple diseases", "Huntington disease", "Oxidative phosphorylation", "Thermogenesis"))[["V2"]]

genes <- unique(c(genes))


subset(M.Promoters.annotated, rownames(M.Promoters.annotated) %in% genes)

genes

M.Promoters.annotated[genes,]

M.Promoters.annotated[c("SIRT6", "ACTL6A"),]

heatmap <- heatmap_func(data              = M.Promoters.annotated[genes,annotations_A1$Sample_Name],
                        annotations       = annotations,
                        subsetAnnotation  = c("Study", "ER_status", "Variant_type", "ATM_LOH"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")


ggsave(plot=heatmap, filename="~/ATM_Analysis/svn/analyse/results/FinalMethylation/VisualisationA1withGenesA2.png", width=12, height=7)










k <- 2
clusters <- cutree(hc, k=k)
target <- ifelse(annotations$Sample_Group == "ATM", 2, 1)
names(target) <- annotations$Sample_Name
makeCluster <- function(hc, k)
{
    clusters <- cutree(hc, k=k)
    confusionM <- table(clusters, target)
    return(confusionM)
}

confusionM <- makeCluster(hc,
                          k = k)
maxInit <- MAX <- max(confusionM)
while (maxInit == MAX)
{
    if(k<9){
    confusionM <- makeCluster(hc,
                              k = k + 1)
    k <- k + 1
    print(confusionM)
    MAX <- max(confusionM)
    print(MAX)
    }
}

clusters <- cutree(hc, k=k)    
confusionM <- table(clusters, target)
# select clusters that have more than half curie samples and others
rows <- rownames(confusionM)[confusionM[, 2] > (confusionM[, 2] + confusionM[, 1])/2]
notRows <- rownames(confusionM)[! rownames(confusionM) %in% rows]
# Replace data in clusters matrix
clusterInRows <- clusters[clusters %in% rows]
clusterInRows <- replace(clusterInRows, clusterInRows %in% rows, 2)
clusterNotInRows <- clusters[clusters %in% notRows]
clusterNotInRows <- replace(clusterNotInRows, clusterNotInRows %in% notRows, 1)
clusters <- c(clusterInRows, clusterNotInRows)
# New confusion matrix
confusionM <- table(clusters, target)

# Metrics
sum(diag(confusionM)) / sum(confusionM)
aricode::NMI(clusters, target, variant = "sum")
mclust::adjustedRandIndex(clusters, target)


dend <- as.dendrogram(hc)
dendextend::labels_colors(dend) <- ""
clusters <- clusters[names(dendextend::labels_colors(dend))]
dendextend::labels_colors(dend) <- clusters


plot(dend)





plot(as.dendrogram(hc))
rect.hclust(hc, k=k
)

clustersCol <- replace(clusters, clusters == 2, "red")
clustersCol <- replace(clustersCol, clusters == 1, "blue")

Var = clusters
varCol = gsub("1","red",Var)                        # convert numbers to colours
varCol = gsub("2","blue",varCol)

dendextend::labels_colors(dend) <- ""


# Visualiser le dendogramme
#as.dendrogram(hc) %>% unclass %>% str
#as.dendrogram(hc)[[1]] %>% dendextend::get_nodes_attr("nodePar")



# Replace the first cluster that contains curie samples with the sum and remove others
confusionM[rows[1], ] <- colSums(confusionM[rows, ])
confusionM <- confusionM[rownames(confusionM)!= rows[-1], ]

# Si plus de deux lignes, même opération mais pour les témoins
confusionM[notRows[1], ] <- colSums(confusionM[notRows, ])
confusionM <- confusionM[rownames(confusionM)!= notRows[-1], ]


sum(diag(confusionM)) / sum(confusionM)


rownames(confusionM)[! rownames(confusionM) %in% rows]


clusters <- replace(clusters, clusters %in% rows, rows[1])

table(clusters, target)


sum(diag(confusionM)) / sum(confusionM)
aricode::NMI(clusters, target, variant = "max")

plot(as.dendrogram(hc))
rect.hclust(hc, k=k)






maxRow <- strtoi(as.character(which.max(rowSums(confusionM))))
minRow <- strtoi(as.character(which.min(rowSums(confusionM))))
while (nrow(confusionM) != 2)
{ 
    confusionM[maxRow, ] <- confusionM[maxRow, ] + confusionM[minRow, ]
    confusionM <- confusionM[-minRow, ]
    clusters <- replace(clusters, clusters == minRow, maxRow)

    maxRow <- strtoi(as.character(which.max(rowSums(confusionM))))
    minRow <- strtoi(as.character(which.min(rowSums(confusionM))))

}


sum(diag(confusionM)) / sum(confusionM)

#clusterNoControl <- names(confusionM[,colMax][confusionM[,ifelse(colMax == "1", "2", "1")] > confusionM[,colMax]])
#clusterToAdd <- rownames(confusionM)[! rownames(confusionM) %in% clusterNoControl]
#clusterToAdd <- names(confusionM[,1][confusionM[,1] == max(confusionM[,1])])
#clusters <- replace(clusters, clusters %in% clusterNoControl, as.character(strtoi(clusterToAdd) + 1))

cm <- table(clusters, target)
sum(diag(cm)) / sum(cm)

1-aricode::NMI(clusters, target, variant = "max")

plot(as.dendrogram(hc))
rect.hclust(hc, k=k)






dend <- as.dendrogram(hc)
COLS <- c("turquoise", "orange")
names(COLS) <- unique(annotations$Sample_Group)
dend <- color_labels(dend, col = COLS[annotations[labels(dend)]])

plot(dend)

plot(as.dendrogram(hc))
rect.hclust(hc, k=4)










library(doParallel)
library(ConsensusClusterPlus)

results = ConsensusClusterPlus(as.matrix(M.Promoters.annotated),
                               maxK = 5,
                               reps = 100,
                               pItem = 0.8,
                               pFeature = 1,
                               clusterAlg = "km",
                               distance = "pearson",
                               seed = 1262118388.71279,
                               plot = "png",
                               writeTable = TRUE)

icl = calcICL(results,
              plot="png")

icl[["clusterConsensus"]]


library(NMF)
res <- nmf(Beta.Promoters.annotated, rank=2:6, nrun = 100, .opt='vp1')

compare(res, class = ifelse(annotations$Sample_Group == "ATM", 1, 0))


consensusmap(res)


s <- extractFeatures(res)





test <- assignGroups_func(matrix = M.Promoters.annotated,
                          sampleTable = annotations,
                          distance          = "pearson",
                          method            = "ward",
                          k                 = 2,
                  #seed              = seed,
                          group_names       = paste0("group_",1:2),
                          group_num_column  = "group_num",
                          group_name_column = "group_name",
                          silhouette        = FALSE
                          )

    





rand.index <- function(group1, group2)
{
    x <- abs(sapply(group1, function(x) x - group1))
    x[x > 1] <- 1
    y <- abs(sapply(group2, function(x) x - group2))
    y[y > 1] <- 1
    sg <- sum(abs(x - y))/2
    bc <- choose(dim(x)[1], 2)
    ri <- 1 - sg/bc
    return(ri)
}






library(biomaRt)
ensembl <-  useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", mirror="uswest")
mart <- useDataset(dataset="hsapiens_gene_ensembl", mart=ensembl)

histoneList <- biomaRt::getBM("hgnc_symbol", filter = "go", values = "GO:0016570", mart = ensembl)


library("org.Hs.eg.db")
histoneList <- AnnotationDbi::select(org.Hs.eg.db, keytype = "GOALL", keys = "GO:0016570", columns = "SYMBOL")


test <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype="ENSEMBL"),
                      columns = c("SYMBOL", "GENENAME"), keytype = "ENSEMBL")

test <- test[grep("histone", test$GENENAME), ]$SYMBOL

toPlot <- M.Promoters.annotated[unique(test),]
toPlot <- toPlot[complete.cases(toPlot),]


plot <- heatmap_func(data          = toPlot,
             annotations   = annotations,
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of genes linked to histones Gene Ontology")
             
ggsave("/data/kdi_prod/.kdi/project_workspace_0/833/acl/05.00/svn/analyse/results/FinalMethylation/Histone_heatmap.png", plot, height = 7, width = 14)


# Mean methylation in genes
M.Genes.mean <- data.frame(t(colMeans(M.Genes.annotated)))
rownames(M.Genes.mean) <- "mean"
plotByGene(data = M.Genes.mean,
           Gene = "mean",
           main = "Mean level of methylation in genes")

# Mean methylation in promoters
M.pro.mean <- data.frame(t(colMeans(M.Promoters.annotated)))
rownames(M.pro.mean) <- "mean"
plotByGene(data = M.pro.mean,
           Gene = "mean",
           main = "Mean level of methylation in promoters")

# Mean methylation in ATM gene
plotByGene(data = M.Genes.annotated,
           Gene = "ATM",
           main = "Level of methylation in ATM gene")

# Mean methylation in ATM promoter
plotByGene(data = M.Promoters.annotated,
           Gene = "ATM",
           main = "Level of methylation in ATM promoter")



# test the relationship between ranks of ATM promoter methylation and general methylation
D1 <- data.frame("ATM" = t(M.Promoters.annotated["ATM", ]), "test" = colMeans(M.Promoters.annotated))
D1 <- subset(D1, rownames(D1) %in% annotations[annotations$Sample_Group == "ATM",]$Sample_Name)
cor.test(x = D1[, 1], y = D1[, 2], method = "spearman")


# test the relationship between ranks of ATM promoter methylation and all other promoters
data <- t(M.Promoters.annotated)
data <- subset(data, rownames(data) %in% annotations[annotations$Sample_Group == "ATM",]$Sample_Name)

#N <- sample(1:dim(data)[2], 800, replace = FALSE)
#data1 <- data[,N]

WilcoxTest <- function(df)
{
    plyr::adply(df, 2, function(x) {
        test <- broom::tidy(cor.test(x = data[, "ATM"],
                                        y = x,
                                        method = "spearman",
                                        exact = FALSE)
                            )
        out <- data.frame(#"var1" = colnames(df)[x[1]],
                          #"var2" = colnames(Data[x[2]]),
            "rho"= test$estimate,
            "statistic"= test$statistic,
            "p.value" = sprintf("%.3f", test$p.value),
            "method"= test$method,
            "alternative"= test$alternative
        )
        return(out)    
    })
}

res <- WilcoxTest(data)

# Filter on rho value
genes <- res[abs(res[,2]) > 0.85, ] 
# Filter on p-value
genes <- genes[genes[,4] < 0.05, ]
genes <- genes[,1]
genes <- as.character(genes)

heatmap_func(data          = M.Promoters.annotated[genes, ],
             annotations   = annotations,
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of promoters which methylation is correlated to ATM promoter")
        


# Mean methylation in ATM, BRCA1 and CDK5 genes
p1 <- plotByGene(data = M.Promoters.annotated,
           Gene = "ATM",
           main = "Level of methylation of ATM promoter")
p2 <- plotByGene(data = M.Promoters.annotated,
           Gene = "BRCA1",
           main = "Level of methylation of BRCA1 promoter")
p3 <- plotByGene(data = M.Promoters.annotated,
           Gene = "CDK5",
           main = "Level of methylation of CDK5 promoter")
ggpubr::ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3,
                  common.legend = TRUE, legend = "right")

# Mean methylation level in H2AX, MDC1 and CK2 promoters (players in the activation of ATM)
p1 <- plotByGene(data = M.Promoters.annotated,
           Gene = "H2AX",
           main = "Level of methylation of H2AX promoter")
p2 <- plotByGene(data = M.Promoters.annotated,
           Gene = "MDC1",
           main = "Level of methylation of MDC1 promoter")
p3 <- plotByGene(data = M.Promoters.annotated,
           Gene = "CSNK2A1",
           main = "Level of methylation of CSNK2A1 promoter")
p4 <- plotByGene(data = M.Promoters.annotated,
           Gene = "CSNK2A2",
           main = "Level of methylation of CSNK2A2 promoter")

ggpubr::ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 4,
                  common.legend = TRUE, legend = "right")

# Mean methylation level in MRN complex
p1 <- plotByGene(data = M.Promoters.annotated,
           Gene = "ATM",
           main = "Level of methylation of ATM promoter")
p2 <- plotByGene(data = M.Promoters.annotated,
           Gene = "MRE11",
           main = "Level of methylation of MRE11 promoter")
p3 <- plotByGene(data = M.Promoters.annotated,
           Gene = "RAD50",
           main = "Level of methylation of RAD50 promoter")
p4 <- plotByGene(data = M.Promoters.annotated,
           Gene = "ATMIN",
           main = "Level of methylation of ATMIN promoter")
ggpubr::ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 4,
                  common.legend = TRUE, legend = "right")

# Mean methylation level in MRN complex
p1 <- plotByGene(data = M.Promoters.annotated,
           Gene = "RB1",
           main = "Level of methylation of RB1 promoter")
p2 <- plotByGene(data = M.Promoters.annotated,
           Gene = "CDK2",
           main = "Level of methylation of CDK2 promoter")
p3 <- plotByGene(data = M.Promoters.annotated,
           Gene = "TP53",
           main = "Level of methylation of TP53 promoter")
p4 <- plotByGene(data = M.Promoters.annotated,
           Gene = "E2F1",
           main = "Level of methylation of E2F1 promoter")
ggpubr::ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 4,
                  common.legend = TRUE, legend = "right")





heatTotal <- heatmap_func(data          = M.Promoters.annotated[, annotations_A1$Sample_Name],
                          annotations   = subset(annotations, annotations$Sample_Name %in% annotations_A1$Sample_Name),
                          fontsize_row  = 8,
                          fontsize_col  = 10,
                          main = "Heatmap of global methylation of promoters")

heatTotal <- heatmap_func(data          = M.Promoters.annotated,
                          annotations   = annotations,
                          fontsize_row  = 8,
                          fontsize_col  = 10,
                          main = "Heatmap of global methylation of promoters")

heatmap_func(data          = M.Promoters.annotated[subset(pro.A1, pro.A1$logFC > 0)$Gene, ],
             annotations   = annotations,
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of genes found with a promoter Differentially methylated")

heatmap_func(data          = M.Promoters.annotated[pro.A1$Gene, ],
             annotations   = annotations,
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of genes found with a promoter Differentially methylated")






genesToPlot <- c("RB1", "CDK2", "TP53", "E2F1", "HDAC1", "H2AX", "CDK4", "TFDP1", "ABL1", "MDM2", "BRCA1", "ATM")

heatmap_func(data          = M.Promoters.annotated[genesToPlot, ],
             annotations   = annotations,
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of promoters which methylation is correlated to ATM promoter")

# Interactors of ATM
genesToPlot <- read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation/ATMandRB1Interactors.csv")[,2]

Matrix <- M.Promoters.annotated[genesToPlot, ][complete.cases(M.Promoters.annotated[genesToPlot, ]), ]
# All samples
heatmap_func(data          = Matrix,
             annotations   = annotations,
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of promoters of genes linked to ATM and RB1")

library("RColorBrewer")

# Pathogenic or likely pathogenic samples
heatmap_func(data          = Matrix[, annotations_A1$Sample_Name],
             annotations   = subset(annotations, annotations$Sample_Name %in% annotations_A1$Sample_Name),
             color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name =
                                                      "RdBu")))(100),
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of promoters of genes linked to ATM and RB1",
             filename        = "~/ATM_Analysis/svn/analyse/results/FinalMethylation/test_scale_heatmap2.png")



rev(c(RColorBrewer::brewer.pal(9, "RdBu")[0], RColorBrewer::brewer.pal(9, "RdBu")[5], RColorBrewer::brewer.pal(9, "RdBu")[9]))



heatmap_func(data          = M.Promoters.annotated[c("AKT2", "MTOR", "TP53"), ],
             annotations   = annotations,
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of promoters of genes linked to ATM and RB1")


#"H4C13", "MIR3187", "GAS7", "MIR3193", "ARHGAP28-AS1", "SCGB3A1", "H3C11", "ANKS6", "ARHGAP28", "ADAMTS12"






AggregatePathwaysPCA <- function(EnrichmentObject,
                              ListOfPathways,
                              type     = "all",
                              freq     = NULL,
                              k        = 2,
                              samples  = NULL # By default all samples are used
)
{
    # get informations for the name of the file
    Feature <- strsplit(EnrichmentObject, "\\.")[[1]][2]
    Analysis <- strsplit(EnrichmentObject, "\\.")[[1]][1]
    name <- paste0(Analysis, ".", Feature)

    # Evaluate the object (from string to object)
    EnrichmentObject <- eval(parse(text = EnrichmentObject))
    # Get the gene names for each pathways indicated as input
    List <- ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description %in% ListOfPathways,],
                               column = "core_enrichment")
    # Handle if you want to filter genes to have genes present in a certain percentage of pathways
    if (type == "intercept")
    {
        Freq <- table(List$core_enrichment)
        if (! is.null(freq))
        {
            Genes <- names(Freq[Freq >= length(ListOfPathways)*freq])
        } else {
            message("Be CAREFULL, freq parameter is not defined and all genes of all found pathways will be used")
            Genes <- List$core_enrichment
        }
        if (length(Genes) == 0)
        {
            stop("Your criteria is too much stringent, no genes are retained and no Heatmap will be created")
        }
    } else {
        Genes <- List$core_enrichment
    }
    # Suppress duplicated genes
    Genes <- unique(Genes)
    #print(Genes)
    print(length(Genes))
    if (dim(List)[1] != 0){
        # Handle if the filename is too long
        fileName <- paste0(unique(List$Description), collapse = "\t")
        fileName <- stringr::str_replace_all(fileName, " - ", "AND")
        fileName <- stringr::str_replace_all(fileName, "-", "_")
        fileName <- stringr::str_replace_all(fileName, " ", "_")
        if (nchar(fileName) > 20)
        {
            fileName <- paste0(substr(fileName, 1, 60))
        }
        # Which Feature is used ==> choose the correct dataset
        if (Feature == "Pro")
        {
            dataset <- M.Promoters.annotated
        }
        if (Feature == "Genes")
        {
            dataset <- M.Genes.annotated
        }
        # handle if you want only a subset of samples
        if (! is.null(samples))
            {
                dataset = dataset[, samples]
                annotations = subset(annotations, annotations$Sample_Name %in% samples)
            }
        # Create and save the heatmap
        main = paste0(unique(List$Description), collapse = " & ")
        PCA <- pca_func(matrix           = dataset,
                        sampleTable      = annotations,
                        sampleLabel      = "Sample_Name",
                        sampleLabelSize  = 4,
                        pointSize        = 5, 
                        highlightedVar   = "Sample_Group",
                        title = main)
        
        ggsave(plot = PCA,
               width = 14,
               height = 7,
               filename = file.path(OutFiles,
                                    paste0(name,
                                           "_PCA_",
                                           ifelse(!is.null(freq), freq, ""),
                                           "_",
                                           ifelse(!is.null(samples),
                                                  stringr::str_split(
                                                               stringr::str_split(
                                                                            deparse(substitute(annotations_A1$Sample_Name)), "_")[[1]][2], "\\$")[[1]][1]
                                                , ""),
                                           fileName,
                                           ".png")
                                    )
               )
    }
    rm(PCA)
}  





############################# Promoters
##### A1
# Cluster 1 with only A1 samples
ListOfPathways <- c("Thermogenesis", "Diabetic cardiomyopathy", "Alzheimer disease", "Huntington disease", "Proteasome", "Pathways of neurodegeneration - multiple diseases", "Parkinson disease", "Prion disease", "Chemical carcinogenesis - reactive oxygen species", "Amyotrophic lateral sclerosis", "Oxidative phosphorylation", "Non-alcoholic fatty liver disease")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name) 

# Cluster 1
ListOfPathways <- c("Thermogenesis", "Diabetic cardiomyopathy", "Alzheimer disease", "Huntington disease", "Proteasome", "Pathways of neurodegeneration - multiple diseases", "Parkinson disease", "Prion disease", "Chemical carcinogenesis - reactive oxygen species", "Amyotrophic lateral sclerosis", "Oxidative phosphorylation", "Non-alcoholic fatty liver disease")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways) 

# Cluster2 with only A1 samples
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name)

# Cluster2
ListOfPathways <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways)

# Cluster 3 with only A1 samples
ListOfPathways <- c("Human papillomavirus infection", "Breast cancer", "Gastric cancer", "Hepatocellular carcinoma", "Pathways in cancer", "EGFR tyrosine kinase inhibitor resistance", "PI3K-Akt signaling pathway", "Glioma", "Proteoglycans in cancer", "Human cytomegalovirus infection")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name)

# Cluster 3
ListOfPathways <- c("Human papillomavirus infection", "Breast cancer", "Gastric cancer", "Hepatocellular carcinoma", "Pathways in cancer", "EGFR tyrosine kinase inhibitor resistance", "PI3K-Akt signaling pathway", "Glioma", "Proteoglycans in cancer", "Human cytomegalovirus infection")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways)


# Cluster 4 with only A1 samples
ListOfPathways <- c("Fanconi anemia pathway", "Homologous recombination")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name)

# Cluster 4
ListOfPathways <- c("Fanconi anemia pathway", "Homologous recombination")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways)


# Cluster 5 with only A1 samples
ListOfPathways <- c("Drug metabolism - cytochrome P450", "Steroid hormone biosynthesis")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways,
                  samples = annotations_A1$Sample_Name)

# Cluster 5
ListOfPathways <- c("Drug metabolism - cytochrome P450", "Steroid hormone biosynthesis")
AggregatePathwaysPCA("A1.Pro.KeggGse", ListOfPathways)









pca_func(matrix           = Matrix,
         sampleTable      = annotations,
         sampleLabel      = "Sample_Name",
         sampleLabelSize  = 4,
         title            = "PCA of genes found within regions with copy number alterations",
         highlightedVar   = "Sample_Group")


tsne_func(matrix           = Matrix,
          sampleTable      = annotations,
          sampleLabel      = "Sample_Name",
          perplexity       = 17,
          fontsize         = 10,
          sampleLabelSize  = 4,
         #sampleLabelSize  = 4,
         #title            = "PCA of genes found within regions with copy number alterations",
          highlightedVar   = "Sample_Group")



MostVariables <- variable_feature_selection(data = M.Promoters.annotated,
                                            thres.num = 5000)

heatmap_func(data          = M.Promoters.annotated,
             AnnotationRB1 = TRUE,
             annotations   = annotations,
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of genes found within regions with copy number alterations")


pca_func(matrix           = M.Promoters.annotated[MostVariables,],
         sampleTable      = annotations,
         sampleLabel      = "Sample_Name",
         sampleLabelSize  = 4,
         title            = "PCA of genes found within regions with copy number alterations",
         highlightedVar   = "Sample_Group")




# Genes with DM promoters found within regions with CN variations 

write.csv2(sapply(M.Promoters.annotated[c("TUSC8", "RNU6-75P", "COX5BP6", "NFATC3", "APOL6"), ], gsub, pattern = ",", replacement= ".")
, file = "~/ATM_Analysis/svn/analyse/results/FinalMethylation/genesInterest.csv")


## Samples with pathogenic or ikely pathogenic variants
heatmap_func(data          = M.Promoters.annotated[c("TUSC8", "RNU6-75P", "COX5BP6", "NFATC3", "APOL6"), annotations_A1$Sample_Name],
             otherAnnotations = TRUE,
             annotations   = subset(annotations, annotations$Sample_Name %in% annotations_A1$Sample_Name),
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of genes found within regions with copy number alterations")


pca_func(matrix           = M.Promoters.annotated[c("TUSC8", "RNU6-75P", "COX5BP6", "NFATC3", "APOL6"), annotations_A1$Sample_Name],
         sampleTable      = subset(annotations, annotations$Sample_Name %in% annotations_A1$Sample_Name),
         sampleLabel      = "Sample_Name",
         sampleLabelSize  = 4,
         title            = "PCA of genes found within regions with copy number alterations",
         highlightedVar   = "Sample_Group")


tsne_func(matrix           = M.Promoters.annotated[c("TUSC8", "RNU6-75P", "COX5BP6", "NFATC3", "APOL6"), annotations_A1$Sample_Name],
          sampleTable      = subset(annotations, annotations$Sample_Name %in% annotations_A1$Sample_Name),
          sampleLabel      = "Sample_Name",
          perplexity       = 17,
          fontsize         = 10,
          sampleLabelSize  = 4,
         #sampleLabelSize  = 4,
         #title            = "PCA of genes found within regions with copy number alterations",
          highlightedVar   = "Sample_Group")


## All samples
heatmap_func(data          = M.Promoters.annotated[c("TUSC8", "RNU6-75P", "COX5BP6", "NFATC3", "APOL6"), ],
             otherAnnotations = TRUE,
             annotations   = annotations,
             fontsize_row  = 8,
             fontsize_col  = 10,
             main = "Heatmap of genes found within regions with copy number alterations")

pca_func(matrix           = M.Promoters.annotated[c("TUSC8", "RNU6-75P", "COX5BP6", "NFATC3", "APOL6"), ],
         sampleTable      = annotations,
         sampleLabel      = "Sample_Name",
         sampleLabelSize  = 4,
         title            = "PCA of genes found within regions with copy number alterations",
         highlightedVar   = "Sample_Group")


tsne_func(matrix           = M.Promoters.annotated[c("TUSC8", "RNU6-75P", "COX5BP6", "NFATC3", "APOL6"), ],
          sampleTable      = annotations,
          sampleLabel      = "Sample_Name",
          perplexity       = 10,
          fontsize         = 10,
          sampleLabelSize  = 4,
         #sampleLabelSize  = 4,
         #title            = "PCA of genes found within regions with copy number alterations",
          highlightedVar   = "Sample_Group")







# Interactors of ATM (found via Cytoscape, query of 100 genes, score of 90%)
genesToPlot <- read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation/ATMandRB1Interactors.csv")[,2]
Matrix <- M.Promoters.annotated[genesToPlot, ][complete.cases(M.Promoters.annotated[genesToPlot, ]), ]

InteractorsATM <- heatmap_func(data          = Matrix,
                               annotations   = annotations,
                               #color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                               fontsize_row  = 8,
                               fontsize_col  = 10,
                               main = "Heatmap of promoters of genes linked to ATM",
                               filename = "~/ATM_Analysis/svn/analyse/results/FinalMethylation/Heatmap_ATM_Interactors_AllSamples.png")

InteractorsATM <- heatmap_func(data          = Matrix[, annotations_A1$Sample_Name],
                               annotations   = subset(annotations, annotations$Sample_Name %in% annotations_A1$Sample_Name),
                               fontsize_row  = 8,
                               fontsize_col  = 10,
                               main = "Heatmap of promoters of genes linked to ATM",
                               filename = "~/ATM_Analysis/svn/analyse/results/FinalMethylation/Heatmap_ATM_interactors_A1Samples.png")









############################# Promoters
##### A1
# Cluster 1 with only A1 samples
ListOfPathwaysGB <- c("Thermogenesis", "Diabetic cardiomyopathy", "Alzheimer disease", "Huntington disease", "Proteasome", "Pathways of neurodegeneration - multiple diseases", "Parkinson disease", "Prion disease", "Chemical carcinogenesis - reactive oxygen species", "Amyotrophic lateral sclerosis", "Oxidative phosphorylation", "Non-alcoholic fatty liver disease")
# Evaluate the object (from string to object)
EnrichmentObject <- eval(parse(text = "A1.Pro.KeggGse"))
# Get the gene names for each pathways indicated as input
PathwaysGroupB <- unique(ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description %in% ListOfPathwaysGB,],
                                  column = "core_enrichment")$core_enrichment)

# Cluster2 with only A1 samples
ListOfPathwaysGA <- c("Systemic lupus erythematosus", "Alcoholism", "Neutrophil extracellular trap formation")
# Evaluate the object (from string to object)
EnrichmentObject <- eval(parse(text = "A1.Pro.KeggGse"))
# Get the gene names for each pathways indicated as input
PathwaysGroupA <- unique(ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description %in% ListOfPathwaysGA,],
                                  column = "core_enrichment")$core_enrichment)

# Cluster 3
ListOfPathwaysGC <- c("Human papillomavirus infection", "Breast cancer", "Gastric cancer", "Hepatocellular carcinoma", "Pathways in cancer", "EGFR tyrosine kinase inhibitor resistance", "PI3K-Akt signaling pathway", "Glioma", "Proteoglycans in cancer", "Human cytomegalovirus infection")
# Evaluate the object (from string to object)
EnrichmentObject <- eval(parse(text = "A1.Pro.KeggGse"))
# Get the gene names for each pathways indicated as input
PathwaysGroupC <- unique(ParseEnrichResults(enrichResultObject = EnrichmentObject[EnrichmentObject$Description %in% ListOfPathwaysGC,],
                                  column = "core_enrichment")$core_enrichment)




AllGenes <- unique(c(PathwaysGroupB, PathwaysGroupA, PathwaysGroupC, genesToPlot))
AllGenes <- AllGenes[order(AllGenes)]
GenesPathways <- data.frame(row.names = AllGenes,
                            "PathwaysGroupB" = AllGenes %in% PathwaysGroupB,
                            "PathwaysGroupA" = AllGenes %in% PathwaysGroupA,
                            "PathwaysGroupC" = AllGenes %in% PathwaysGroupC,
                            "ATM interactors" = AllGenes %in% genesToPlot)

GenesPathways[GenesPathways == TRUE] <- "Included"
GenesPathways[GenesPathways == FALSE] <- "Not included"


knitr::kable(GenesPathways, "html", booktabs = TRUE, align = c('l', rep('c', 3)))  %>%
    kableExtra::kable_styling(latex_options = "striped", bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F) %>%
        kableExtra::save_kable("~/ATM_Analysis/svn/analyse/results/FinalMethylation/Interpretation/Kegg/Gse/PromotersPathways.html")



queryGene <- function(enrichResultObject,
                      comparison = NULL,
                      gene
                      )
    {
        EnrichmentObject <- eval(parse(text = "A1.Pro.KeggGse"))
        if (!is.null(comparison))
            {
                EnrichmentObject <- EnrichmentObject[EnrichmentObject$comparison == comparison,]
            }
        rank <- grepl(gene, strsplit(EnrichmentObject[, "core_enrichment"], "/"))
        return(A1.Pro.KeggGse$Description[rank])
    }


queryGene("A1.Pro.KeggGse", "A1", "ATM")

ListOfPathwaysGB[ListOfPathwaysGB %in% queryGene("A1.Pro.KeggGse", "A1", "ABL1")]

ListOfPathwaysGB[ListOfPathwaysGB %in% queryGene("A1.Pro.KeggGse", "A1", "TP53")]


ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "ATM")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "ATR")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "BRCA1")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "BRCA2")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "CDK2")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "CDK4")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "MLH1")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "RAD51")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "STK11")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "TP53")]

ListOfPathwaysGC[ListOfPathwaysGC %in% queryGene("A1.Pro.KeggGse", "A1", "TRAF6")]






CA.A1 <- c("T0072", 
           "T0076", 
           "T0077R", 
           "T0078", 
           "T0123", 
           "T0181", 
           "T0247", 
           "T0252")

CA.All <- c("T0076",
            "T0077R",
            "T0078",
            "T0100",
            "T0106",
            "T0109",
            "T0120",
            "T0123",
            "T0181",
            "T0220",
            "T0232",
            "T0247",
            "T0248",
            "T0249",
            "T0252")

CB.A1 <- c("T0072",
            "T0076",
            "T0078",
            "T0091",
            "T0123",
            "T0181",
            "T0220",
            "T0247",
            "T0248")

CB.All <- c("T0072",
            "T0076",
            "T0077R",
            "T0078",
            "T0091",
            "T0100",
            "T0106",
            "T0109",
            "T0120",
            "T0123",
            "T0141",
            "T0181",
            "T0220",
            "T0232",
            "T0247",
            "T0248",
            "T0249")

CC.A1 <- c("T0009",
           "T0072",
           "T0076",
           "T0077R",
           "T0078",
           "T0091",
           "T0123",
           "T0141",
           "T0181",
           "T0220",
           "T0247",
           "T0248",
           "T0249")

CC.All <- c("T0072",
            "T0076",
            "T0077R",
            "T0078",
            "T0091",
            "T0100",
            "T0106",
            "T0109",
            "T0120",
            "T0123",
            "T0141",
            "T0181",
            "T0220",
            "T0232",
            "T0247",
            "T0248")

samples <- annotations$Sample_Name[startsWith(annotations$Sample_Name, "T")]
comp <- data.frame(row.names = samples,
                   "Cluster Group A - A1" = samples %in% CA.A1,
                   "Cluster Group A - All" = samples %in% CA.All,
                   "Cluster Group B - A1" = samples %in% CB.A1,
                   "Cluster Group B - All" = samples %in% CB.All,
                   "Cluster Group C - A1" = samples %in% CC.A1,
                   "Cluster Group C - All" = samples %in% CC.All
                   )

comp[comp == TRUE] <- "Included"
comp[comp == FALSE] <- "Not included"

removedSamples <- rownames(initial_annotations)[!rownames(initial_annotations) %in% rownames(annotations)]
notUsedSamples <- removedSamples[!startsWith(removedSamples, "T")]
removedSamples <- removedSamples[startsWith(removedSamples, "T")]

plot_comp <- PlotAnnotation(annotations     = comp,
                            brewer          = TRUE,
                            sameColors      = TRUE,
                            transpose       = TRUE
                            )

ggsave(plot = plot_comp,
       filename = file.path(OutFiles, "InclusionInAnalyses.png"),
       width = 60, height = 1/3*60,
       units = "cm")












eeee



























plotByGene(data = M.Genes.annotated,
           Gene = "ATM",
           main = "Level of methylation in ATM gene")

plotByGene(data = M.Promoters.annotated,
           Gene = "ATM",
           main = "Level of methylation in ATM promoter")


plotByGene(data = M.Genes.annotated,
           Gene = "BRCA1",
           main = "Level of methylation in BRCA1 gene")

plotByGene(data = M.Promoters.annotated,
           Gene = "BRCA1",
           main = "Level of methylation in BRCA1 promoter")


plotByGene(data = M.Genes.annotated,
           Gene = "BRCA2",
           main = "Level of methylation in BRCA2 gene")

plotByGene(data = M.Promoters.annotated,
           Gene = "BRCA2",
           main = "Level of methylation in BRCA2 promoter")





















############ If for all analyses


# getAllGenesForOneEnrichedPathway: function to recover all genes found in different analysis in a particular pathway
# @DataFrames: dataframe in which to look for
# @column: in which column are indicated the genes
getAllGenesForOneEnrichedPathway <- function(DataFrames = c("A1.Pro.KeggGse",
                                                            "A2.Pro.KeggGse",
                                                            "A3.Pro.KeggGse",
                                                            "A4.Pro.KeggGse"),
                                             column = "core_enrichment")
{
    dict = c()
    for (N in DataFrames)
    {
        for (i in eval(parse(text=N))$Description)
        {
            dataFrame = eval(parse(text=N))
            listGenes <- ParseEnrichResults(enrichResultObject = dataFrame[dataFrame$Description == i[1],],
                                            column = column)
            dict[[i]] <- c(dict[[i]], listGenes)
            dict[[i]] <- dict[[i]][!duplicated(dict[[i]])]
            dict[[i]] <- dict[[i]][order(dict[[i]])]
        }
    }
    return(dict)
}


getAllGenes <- function(objectPathways,
                        pathways)
{
    list = c()
    for (p in pathways)
    {
        list <- c(list, objectPathways[[p]])
    }
    return(list)
}




#all.ProGSEA <- getAllGenesForOneEnrichedPathway()

#ListOfPathways <- c("Transcriptional misregulation in cancer", "Basal transcription factors", "Fanconi anemia pathway")
#List <- getAllGenes(all.ProGSEA, ListOfPathways) 


#heatmap_func(data          = M.Promoters.annotated[List,],
#                                    annotations   = annotations,
#                                    fontsize_row  = 8,
#                                    fontsize_col  = 10,
#                                    main = paste(ListOfPathways, collapse = " & "),
#                                    filename = "~/ATM_Analysis/svn/analyse/results/methylation/Heatmap_Pro_Transcriptional_miregulation_Basal_transcription_Fanconi.png")


#k <- 2
#umap.plot <- umap_func(matrix              = M.Promoters.annotated[List,],
#                       sampleTable         = annotations,
#                       highlightedVar      = "Sample_Group",
#                       highlightedVarName  = "Sample_Group",
#                       colors              = c("orange", "Red", "Yellow", "Green", "Blue"),
#                       sampleLabel         = "Sample_Name",
#                       metric              = "euclidean",
#                       title               = paste(ListOfPathways, collapse = " & "),
#                       fontsize            = 14,
#                       assignGroup         = TRUE,
#                       )
#ggsave(filename = "~/ATM_Analysis/svn/analyse/results/methylation/Umap_Pro_Transcriptional_miregulation_Basal_transcription_Fanconi.png", plot = umap.plot, width = 12, height = 7)











##k <- 2
#umap_func(matrix              = M.Promoters.annotated[List,],
#          sampleTable         = annotations,
#          highlightedVar      = "Sample_Group",
#          highlightedVarName  = "Sample_Group",
#          colors              = c("orange", "Red", "Yellow", "Green", "Blue"),
#          sampleLabel         = "Sample_Name",
#          metric              = "euclidean",
#          title               = "UMAP",
#          fontsize            = 14,
#          assignGroup         = TRUE,
#          )


#assignGroups_func(as.matrix(M.Promoters.annotated[List, ]),
#                  sampleTable = annotations,
#                  distance          = "euclidean",
#                  method            = "ward.D2",
#                  k                 = 2,
#                  group_names       = paste0("group_",1:k),
#                  group_num_column  = "Sample_Group",
#                  group_name_column = c("ATM", "Non.ATM"),
#                  silhouette        = FALSE
#                  )            


head(M.Promoters.annotated)


test <- heatmap_func(data = M.Promoters.annotated[1:10,],
             annotations       = annotations,
             subsetAnnotation  = c("Study", "ER_status", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )


Histones <- rownames(M.Promoters.annotated)[grep("^H[2|4]", rownames(M.Promoters.annotated))]


test <- heatmap_func(data = M.Promoters.annotated[Histones,],
             annotations       = annotations,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )





test <- heatmap_func(data = M.Promoters.annotated[c("CD86", "CD40", "PROM1"),],
             annotations       = annotations,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )





vec <- c(2,3,4)

names(test[2])

names(test[[names(test[2])]])

for (i in seq_along(names(test)))
    {
        names(test[[names(test[i])]])[names(test[[names(test[i])]]) == " "] <- NA
    }


ComplexHeatmap::decorate_annotations("points")


test

col_fun = circlize::colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
ha = ComplexHeatmap::HeatmapAnnotation(
    foo = 1:10, 
    bar = sample(letters[1:3], 10, replace = TRUE),
    col = list(foo = col_fun,
               bar = c("a" = "red", "b" = "green", "c" = "blue")
    ))
ComplexHeatmap::draw(ha)






set.seed(123)
ha1 = ComplexHeatmap::HeatmapAnnotation(df = data.frame(type = rep(letters[1:2], 5),
                                                        type2 = rep(letters[1:2], 5)))
ComplexHeatmap::Heatmap(matrix(rnorm(100), 10), name = "mat", km = 2,
            top_annotation = ha1)

ComplexHeatmap::decorate_annotation("type",
{
    grid::grid.points(x = unit(c(0.25, 0.45, 0.6, 0.8), "npc"),
                      y = rep(unit(c(0.5), "npc"), 4),
#                      r = 0.25,
                      pch = 8
                      )
})

ComplexHeatmap::decorate_annotation("point", {
        grid::grid.rect(gp = grid::gpar(fill = "#FF000080"))
        }, slice = 2)

test





M.Promoters.annotated[test, ]


map <- heatmap_func(data = M.Promoters.annotated[comm, annotations_A1$Sample_Name],
             annotations       = annotations[annotations$Sample_Name %in% annotations_A1$Sample_Name,],
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )



map <- heatmap_func(data = M.Promoters.annotated[comm, annotations_A2$Sample_Name],
             annotations       = annotations[annotations$Sample_Name %in% annotations_A2$Sample_Name,],
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )




test1 <- c("MIR3187",
"EDEM2",
"ING4",
"ZBTB12BP",
"CWC27",
"MED20",
"ZRANB2",
"INO80B",
"INO80B-WBP1",
"WASF1",
"LIN7C",
"UFL1-AS1",
"ZRANB2-DT",
"DCAKD",
"WDR75",
"TAF13",
"NUP205",
"BDNF-AS",
"ISCU",
"TAX1BP1-AS1",
"RPS2P32",
"HIPK2",
"RIMKLBP2",
"VEGFB",
"ST13P5",
"GRIA2",
"GPR27",
"CASTOR3",
"SERPINB5",
"NUTF2P4",
"NACC2",
"COX5BP6",
"RNU7-84P",
"GPC1",
"FGFR3",
"MIR4516",
"RASGEF1C",
"ACTG1P23",
"CASKIN1",
"NCS1",
"DNAJC7")



comm <- c("MIR3187", "ZBTB12BP",
"ING4",
"EDEM2",
"NUP205",
"CWC27",
"ZRANB2",
"UFL1-AS1",
"WASF1",
"INO80B",
"INO80B-WBP1",
"MED20",
"LIN7C",
"ZRANB2-DT",
"RPS2P32",
"DCAKD",
"GPATCH8",
"BDNF-AS",
"ISCU",
"RPUSD1",
"RNU6-75P",
"ANGPTL3",
"NUTF2P4",
"NAT8L",
"GPR27",
"HIPK2",
"VEGFB",
"GPC1",
"COX5BP6",
"GRIA2",
"RNU7-84P",
"NACC2",
"RASGEF1C",
"CASKIN1",
"ACTG1P23",
"FGFR3",
"NCS1",
"MIR4516",
"DNAJC7")



# signal-to-noise (S2N) metric ratio
# Doc: https://haifengl.github.io/api/java/smile/feature/SignalNoiseRatio.html
signalToNoise <- function(data = M.Promoters.annotated,
                          list1,
                          list2,
                          numberOfMarker = 0.2
                          )
{
    data <- data[c(list1, list2)]
    data$mean.l1 <- rowMeans(data[list1]) # Moyenne l1
    data$mean.l2 <- rowMeans(data[list2]) # Moyenne l2
    data$var.l1 <- sqrt(apply(data[list1], 1, var)) # Standard deviation l1
    data$var.l2 <- sqrt(apply(data[list2], 1, var)) # Standard deviation l2

    # Ensure that variance is at least 20% of mean
    #thresholdSD <- function(mean, SD, percent)
    #{
    #    for (i in 1:length(mean)){
    #    returnValue = SD[i]
    #    absMean = abs(mean[i])
    #    minStdev = percent * absMean
    #    if(minStdev > SD[i])
    #    {
    #        returnValue = minStdev
    #    }
    #    if(returnValue < percent)
    #    {
    #        returnValue = percent
    #    }
    #    return(returnValue)
    #    }
    #}
    #data$SD.l1 <- mapply(thresholdSD, data$mean.l1, data$var.l1, percent)
    #data$SD.l2 <- mapply(thresholdSD, data$mean.l2, data$var.l2, percent)
    #data$SignalNoiseRatio <- (data$mean.l1 - data$mean.l2) / (data$SD.l1 + data$SD.l2)
    data$Signal <- (data$mean.l1 - data$mean.l2) / (data$var.l1 + data$var.l2)
    data <- data[order(data['Signal'], decreasing=T),]
    return(rbind(head(data, n=as.integer(numberOfMarker/2)),
                 tail(data, n=as.integer(numberOfMarker/2))))
}

test <- signalToNoise(data = M.Promoters.annotated,
                      list1 = annotations[annotations$Sample_Group == "ATM", ]$Sample_Name,
                      list2 = annotations[annotations$Sample_Group != "ATM", ]$Sample_Name,
                      numberOfMarker = 40
              )




map <- heatmap_func(data = M.Promoters.annotated[rownames(test), annotations_A1$Sample_Name],
             annotations       = annotations[annotations$Sample_Name %in% annotations_A1$Sample_Name,],
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Variant_type", "ATM_LOH", "Grade"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )


testA2 <- signalToNoise(data = M.Promoters.annotated,
                      list1 = annotations_A2[annotations_A2$Sample_Group == "ATM", ]$Sample_Name,
                      list2 = annotations_A2[annotations_A2$Sample_Group != "ATM", ]$Sample_Name,
                      percent = 0.2
              )


map <- heatmap_func(data = M.Promoters.annotated[rownames(testA2)[1:40], annotations_A2$Sample_Name],
             annotations       = annotations[annotations$Sample_Name %in% annotations_A2$Sample_Name,],
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )


sel <- rownames(testA2)[1:40][rownames(testA2)[1:40] %in% rownames(test)[1:40]]

map <- heatmap_func(data = M.Promoters.annotated[sel, annotations_A2$Sample_Name],
             annotations       = annotations[annotations$Sample_Name %in% annotations_A2$Sample_Name,],
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )


testA2




d <- M.Promoters.annotated[c("MECOM", "ETS1", "IRF7", "LYN", "PDXK", "PTGS2", "RUNX3", "VIM", "DAB2IP", "HSD17B4", "PER1"),]


d <- M.Genes.annotated[c("MECOM", "ETS1", "IRF7", "LYN", "PDXK", "PTGS2", "RUNX3", "VIM", "DAB2IP", "HSD17B4", "PER1", "ACADL", "ADAMTSL1", "ARFGAP3", "B3GAT1", "CDCA7", "FAM78A", "FAM89A", "RNF145", "SLFN11", "HAAO", "HEY2", "HOXB9", "ITGA11", "PROX1", "PSAT1", "RECK", "SMOC1", "SND1", "TNFSF9", "ADHFE1", "DYNLRB2", "HSD17B8", "PISD", "WNK4", "PDXK"),]


map <- heatmap_func(data = d,
             annotations       = annotations,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )






test <- umap_func(M.Promoters.annotated[annotations_A1$Sample_Name],
          sampleTable = annotations_A1,
          highlightedVar      = "Sample_Group",
          highlightedVarName  = "Sample_Group",
          colors              = c("orange", "Red", "Yellow", "Green", "Blue"),
          metric              = "euclidean",
          title               = "UMAP")



OutFiles <- "~/ATM_Analysis/svn/analyse/results/FinalMethylation/differentialAnalysis/GenomeWide"


select <- c("T0141", "T0077R", "T0078", "T0072", "T0091", "T0247", "T0181", "T0123", "T0248")

    
annotations_select <- annotations[annotations$Sample_Name %in% select | annotations$Sample_Group == "Non.ATM", ]

DifferentialAnalysis(data            = M.Promoters.annotated[annotations_select$Sample_Name],
                     annotations     = annotations_select,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "TRUE_FALSE_ATM_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )



TF <- read.csv(file.path(OutFiles, "TRUE_FALSE_ATM_Promoters_DMPs.csv"), sep=",", row.names = 1)

geneSel <- c(rownames(head(TF[TF$DifferentiallyMethylated == "TRUE",], n=20)),
             rownames(tail(TF[TF$DifferentiallyMethylated == "TRUE",], n=20)))

rownames(test)



map <- heatmap_func(data = M.Promoters.annotated[geneSel, annotations_select$Sample_Name],
             annotations       = annotations_select,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Grade", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )



TF <- read.csv(file.path(OutFiles, "TRUE_FALSE_ATM_Promoters_DMPs.csv"), sep=",")
names(TF)[1] <- "Gene"
TF <- TF[c("Gene", "logFC", "P.Value", "adj.P.Val")]
test <- getKeggGSE(TF, comparison = "A1", map = FALSE)

print("A1 done")


test <- getKeggGSE(data.genes.A1, comparison = "A1", map = TRUE)



library(BiocParallel)
register(SerialParam())
getKeggGSE <- function(geneList,
                       comparison = NULL,
                       map        = FALSE,
                       seed       = NULL)
{
    # Set seed
    #set.seed(seed)
    feature <- gsub(".*\\.(.*)\\..*", "\\1", deparse(substitute(geneList)))
    
    # Create the dataframe needed for GSEA analysis
    corr <- bitr(geneList$Gene,
                 fromType  = "SYMBOL",
                 toType    = "ENTREZID",
                 OrgDb     = "org.Hs.eg.db")
    d <- merge(corr, geneList, by.x = "SYMBOL", by.y = "Gene")[c("ENTREZID", "logFC")]
    #print(d)
    print("passed")
    geneList = d[,2]
    names(geneList) = as.character(d[,1])
    geneList = sort(geneList, decreasing = TRUE)
        print("passed2")
    # Perform GSEA analysis
    try(kk2 <- gseKEGG(geneList       = geneList,
                       organism       = 'hsa',
                       keyType        = "ncbi-geneid",
                       minGSSize      = 2,
                       maxGSSize      = 500,
                       pvalueCutoff   = 0.05,
                       pAdjustMethod  = "BH",
                       verbose        = FALSE,
                       seed           = TRUE))
    print(kk2)
    print("passed3")
    # Return the dataframe
    if (exists("kk2"))
    {
        if (nrow(data.frame(kk2)) != 0)
        {
            print("here")
            Gse_Kegg <- setReadable(kk2,
                                    OrgDb   = org.Hs.eg.db,
                                    keyType ="ENTREZID")
            print(Gse_Kegg)
    print("passed4")
            Gse_Kegg <- Gse_Kegg %>%
                arrange(desc(abs(NES)))
    print("passed5")
            #slot(Gse_Kegg, "result")$NES_sign <- ifelse(slot(Gse_Kegg, "result")$NES < 0, 'hypomethylated', "hypermethylated")
            return(Gse_Kegg)
            
        } else {
            return(data.frame(matrix(ncol = 0, nrow = 0)))
        }
    }
}

res <- ParseEnrichResults(test)



map <- heatmap_func(data = M.Promoters.annotated[res[res$Description == "p53 signaling pathway",]$core_enrichment,],
             annotations       = annotations_select,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Grade", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )














heatmap_func(data = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("CSTB",
                                                                                 "RPS6KA1",
                                                                                 "OR52H1",
                                                                                 "SMG1P3",
                                                                                 "PDIA3P2",
                                                                                 "POLR2L",
                                                                                 "HNRNPA1P49",
                                                                                 "UBAP2",
                                                                                 "PRR4",
                                                                                 "MALT1-AS1",
                                                                                 "GPR161",
                                                                                 "SLC7A5P2",
                                                                                 "PRH1",
                                                                                 "RPL36AP30",
                                                                                 "PKN2-AS1",
                                                                                 "DSP-AS1",
                                                                                 "ZNF212",
                                                                                 "TAS2R7",
                                                                                 "SMIM10L1",
                                                                                 "PTDSS2"
                                                                                 ), ],
             annotations       = annotations_A1,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Grade", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )




heatmap_func(data = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c('MALT1-AS1', 'FGD3', 'PEBP1P1', 'PDIA3P2', 'TTLL9', 'VN1R5', 'RASGEF1C', 'SGMS1', 'ANXA3', 'SCAPER', 'NUDT13', 'BAZ2B-AS1', 'NELFB', 'AMPD3', 'DUSP15', 'CYFIP1', 'IGKV1OR2-108', 'MMP28', 'BNIP2', 'SLC6A4', 'SIRT7', 'MIR4734', 'FLT4', 'RPL36AP30', 'ARF4', 'KCNK15', 'BRD1', 'ZNF276', 'GBP2', 'OR1S1', 'PHB1P11', 'MIR3936HG', 'ZDHHC20', 'HBAP1', 'LINC00958', 'RPS10P7', 'CBR1', 'SLC22A5', 'TMC1', 'SF3B4', 'AQP12A', 'MIR218-2', 'RPL5P24', 'BMPR1AP2', 'PTPRVP', 'H3P36', 'NCOA4', 'RNA5SP354', 'KCTD4', 'MIR4530', 'NPEPPSP1', 'ZYG11B', 'PKN2-AS1', 'GNPTAB', 'PTDSS2', 'ATP6V0CP3', 'H3-3A', 'INTS6P1', 'CSTB', 'RNU6-72P'), ],
             annotations       = annotations_A1,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Grade", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )


heatmap_func(data = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c('SCAPER', 'MALT1-AS1', 'PDIA3P2', 'FGD3', 'RPL36AP30', 'RASGEF1C', 'BNIP2', 'ARF4', 'PRR4', 'TAS2R14', 'PRH1', 'FLT4', 'NELFB', 'AMPD3', 'DSP-AS1', 'RPS6KA1', 'IGKV1OR2-108', 'TTLL9', 'MRPL16', 'RIPK1', 'RFX1', 'PEBP1P1', 'PPIC-AS1', 'SMG1P3', 'SLC7A5P2', 'SMIM10L1', 'SNORA14A', 'TNS2-AS1', 'HNRNPA1P49', 'CASTOR3P', 'DMBX1', 'KCTD4', 'TM2D1', 'NCOR2', 'PRELID1', 'BAG3', 'ZNF232-AS1', 'NUP43', 'HK2', 'ABCD1P2', 'EMSLR', 'API5P2', 'RPS10P7', 'RN7SKP114', 'MFHAS1', 'AFDN-DT', 'ZYG11B', 'S100PBP', 'LINC01600', 'EIF3B', 'GPR161', 'H3-3A', 'LGALSL', 'PTDSS2', 'GNPTAB', 'PTPRVP', 'INTS6P1', 'PKN2-AS1', 'CSTB', 'RNU6-72P'), ],
             annotations       = annotations_A1,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Grade", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )





FC.threshold     = 1
p.val.threshold  = 0.05
data.pro.A1$"DifferentiallyMethylated" <- ifelse(data.pro.A1$adj.P.Val <= p.val.threshold & abs(data.pro.A1$logFC) >= FC.threshold, "TRUE", "FALSE")
data.pro.A2$"DifferentiallyMethylated" <- ifelse(data.pro.A2$adj.P.Val <= p.val.threshold & abs(data.pro.A2$logFC) >= FC.threshold, "TRUE", "FALSE")

length(data.pro.A1[data.pro.A1$DifferentiallyMethylated == TRUE, ]$Gene)

heatmap_func(data = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% data.pro.A1[data.pro.A1$DifferentiallyMethylated == TRUE, ]$Gene, ],
             annotations       = annotations,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Grade", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )

pca_func(matrix = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% data.pro.A1[data.pro.A1$DifferentiallyMethylated == TRUE, ]$Gene, annotations_A1$Sample_Name],
         sampleTable = annotations_A1,
         highlightedVar = "Variant_type_aggregated")

A1.probes <- read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation/differentialAnalysis/GenomeWide/A1_Probes_DMPs.csv", row.names=1)

A2.probes <- read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation/differentialAnalysis/GenomeWide/A2_Probes_DMPs.csv", row.names=1)

FC.threshold     = 2
p.val.threshold  = 0.05
A1.probes$"DM" <- ifelse(A1.probes$adj.P.Val <= p.val.threshold & abs(A1.probes$logFC) >= FC.threshold, "TRUE", "FALSE")
A2.probes$"DM" <- ifelse(A2.probes$adj.P.Val <= p.val.threshold & abs(A2.probes$logFC) >= FC.threshold, "TRUE", "FALSE")


M.probes <- read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation/matrices/M_Funnorm_Normalisation.csv", row.names = 1)


#rownames(M.probes) <- M.probes$X
#M.probes <- M.probes[, !names(M.probes) %in% c("X")]


length(rownames(A1.probes[A1.probes$DifferentiallyMethylated == TRUE, ]))

heatmap_func(data = M.probes[rownames(M.probes) %in% rownames(A1.probes[A1.probes$DifferentiallyMethylated == TRUE, ]), annotations_A1$Sample_Name],
             annotations       = annotations_A1,
             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Grade", "Variant_type", "ATM_LOH"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )







DirMatrices = "~/ATM_Analysis/svn/analyse/results/FinalMethylation/matrices"
files <- list.files(path = ".", pattern="_DMPs.csv", recursive = FALSE)


A1.before <- read.csv("~/ATM_Analysis/svn/analyse/results/MethylationWithoutAutralia/FinalMethylation/differentialAnalysis/GenomeWide/A1_Promoters_DMPs.csv", row.names=1)
A1.now <- read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation/differentialAnalysis/GenomeWide/A1_Promoters_DMPs.csv", row.names=1)

A1.before <- A1.before[c("logFC", "adj.P.Val")]
A1.now <- A1.now[c("logFC", "adj.P.Val")]

FC.threshold     = 1
p.val.threshold  = 0.05
A1.before$"DM" <- ifelse(A1.before$adj.P.Val <= p.val.threshold & abs(A1.before$logFC) >= FC.threshold, "TRUE", "FALSE")
A1.now$"DM" <- ifelse(A1.now$adj.P.Val <= p.val.threshold & abs(A1.now$logFC) >= FC.threshold, "TRUE", "FALSE")

dim(A1.before[A1.before$DM == TRUE,])
dim(A1.now[A1.now$DM == TRUE,])


table(rownames(A1.before[A1.before$DM == TRUE,]) %in% rownames(A1.now[A1.now$DM == TRUE,]))

table(rownames(A1.now[A1.now$DM == TRUE,]) %in% rownames(A1.before[A1.before$DM == TRUE,]))

# Only consider the p-value
table(rownames(A1.now[A1.now$adj.P.Val <= p.val.threshold,]) %in% rownames(A1.before[A1.before$adj.P.Val <= p.val.threshold,]))

table(rownames(A1.before[A1.before$adj.P.Val <= p.val.threshold,]) %in% rownames(A1.now[A1.now$adj.P.Val <= p.val.threshold,]))


rownames(A1.now[A1.now$adj.P.Val <= p.val.threshold,])


length(rownames(A1.before) %in% rownames(A1.now))
length(rownames(A1.now) %in% rownames(A1.before))

rownames(A1.before[A1.before$adj.P.Val <= p.val.threshold,]) %in% rownames(A1.now[A1.now$adj.P.Val <= p.val.threshold, ])



plotRanks <- function(a, b, labels.offset=0.1, arrow.len=0.1)
{
    old.par <- par(mar=c(1,1,1,1))
    # Find the length of the vectors
    len.1 <- length(a)
    len.2 <- length(b)
    # Plot two columns of equidistant points
    plot(rep(1, len.1), 1:len.1, pch=20, cex=0.8,
         xlim=c(0, 3), ylim=c(0, max(len.1, len.2)),
         axes=F, xlab="", ylab="") # Remove axes and labels
    points(rep(2, len.2), 1:len.2, pch=20, cex=0.8)

              # Put labels next to each observation
    text(rep(1-labels.offset, len.1), 1:len.1, a)
    text(rep(2+labels.offset, len.2), 1:len.2, b)
    
              # Now we need to map where the elements of a are in b
              # We use the match function for this job
    a.to.b <- match(a, b)

              # Now we can draw arrows from the first column to the second
    arrows(rep(1.02, len.1), 1:len.1, rep(1.98, len.2), a.to.b,
           length=arrow.len, angle=20)
    par(old.par)
}


# Only considering p-value
plotRanks(a=rownames(A1.before[A1.before$adj.P.Val <= p.val.threshold,]),
          b=rownames(A1.now[A1.now$adj.P.Val <= p.val.threshold, ]))


table(rownames(A1.before[A1.before$adj.P.Val <= p.val.threshold,]) %in% rownames(A1.now[A1.now$adj.P.Val <= p.val.threshold, ]))



head(A1.before[order(A1.before$adj.P.Val), ])

head(A1.now[order(A1.now$adj.P.Val), ], n=100)

plotRanks(a=rownames(A1.before[0:100,]),
          b=rownames(A1.now[0:100,]))

plotRanks(a=rownames(A1.before[0:500,]),
          b=rownames(A1.now[A1.before$DM == TRUE,][0:267,]))

plotRanks(a=rownames(A1.before[A1.before$DM == TRUE,]),
          b=rownames(A1.now[A1.now$DM == TRUE,]))







#TCGA <- read.csv("/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation_TCGA/matrices/M_Funnorm_Normalisation.csv", row.names=1)

###load("/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation_TCGA/matrices/M_Funnorm_Normalisation.RData") !!!!!! those are unnormalised data !!!!!!
all_M_Funnorm <- read.csv("/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation_TCGA/matrices/M_Funnorm_Normalisation.csv")
rm(all_M_Funnorm)

TCGA.annotations <- read.csv("/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation_TCGA/matrices/Annotations.csv", row.names=1)
rownames(TCGA.annotations) <- TCGA.annotations$barcode
rownames(TCGA.annotations) <- gsub("-", ".", rownames(TCGA.annotations))


rownames(A1.probes[A1.probes$DM == TRUE,])

TCGA[rownames(A1.probes[A1.probes$DM == TRUE,]),]

heatmap_func(data = TCGA[rownames(A1.probes[A1.probes$DM == TRUE,]),],
             annotations = TCGA.annotations,
             subsetAnnotation  = c("Sample_Group"),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'test pheatmap',
             )


library(pheatmap)


data <- TCGA[rownames(A1.probes[A1.probes$DM == TRUE,]),][complete.cases(TCGA[rownames(A1.probes[A1.probes$DM == TRUE,]),]),]



pheatmap(as.matrix(data), annotation_col = TCGA.annotations["Sample_Group"])

,
        annotations = TCGA.annotations,
        subsetAnnotation  = c("Sample_Group"),
        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
        fontsize_row      = 8,
        fontsize_col      = 10,
        main              = 'test pheatmap',
             




data_Promoters <- read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation_TCGA/matrices/M_Promoters_annotated.csv", row.names=1)





    d <- TCGA.annotations[c("Sample_Group")]

    merge(t(data.frame(TCGA.annotations[c("Sample_Group")])))
    

head(data.frame(TCGA.annotations[c("Sample_Group")]))

head(t(data.frame((data_Promoters["PKN2-AS1",]))))



plotByGene <- function(gene){    
    meth <- merge(x = t(data.frame((data_Promoters[gene,]))),
                  y = data.frame(TCGA.annotations[c("Sample_Group")]),
                  by = "row.names")
    names(meth) <- c("Sample_Name", "Methylation", "Sample_Group")

    meth$Methylation <- as.numeric(meth$Methylation)
    
    plot <- ggplot(meth, aes(x=Sample_Group, y = Methylation)) +
        geom_point(position=position_dodge(0.8))
    return(plot)
}

plotByGene("ATM")



eeeee





pheatmap::pheatmap(mat = data_Promoters[rownames(data_Promoters) %in% c('MALT1-AS1', 'FGD3', 'PEBP1P1', 'PDIA3P2', 'TTLL9', 'VN1R5', 'RASGEF1C', 'SGMS1', 'ANXA3', 'SCAPER', 'NUDT13', 'BAZ2B-AS1', 'NELFB', 'AMPD3', 'DUSP15', 'CYFIP1', 'IGKV1OR2-108', 'MMP28', 'BNIP2', 'SLC6A4', 'SIRT7', 'MIR4734', 'FLT4', 'RPL36AP30', 'ARF4', 'KCNK15', 'BRD1', 'ZNF276', 'GBP2', 'OR1S1', 'PHB1P11', 'MIR3936HG', 'ZDHHC20', 'HBAP1', 'LINC00958', 'RPS10P7', 'CBR1', 'SLC22A5', 'TMC1', 'SF3B4', 'AQP12A', 'MIR218-2', 'RPL5P24', 'BMPR1AP2', 'PTPRVP', 'H3P36', 'NCOA4', 'RNA5SP354', 'KCTD4', 'MIR4530', 'NPEPPSP1', 'ZYG11B', 'PKN2-AS1', 'GNPTAB', 'PTDSS2', 'ATP6V0CP3', 'H3-3A', 'INTS6P1', 'CSTB', 'RNU6-72P'), ],
             annotation_col = TCGA.annotations["Sample_Group"],
             cluster_distance_col = "correlation",
             cluster_distance_row = "correlation")






pheatmap::pheatmap(mat = data_Promoters[rownames(data_Promoters) %in% c('SCAPER', 'MALT1-AS1', 'PDIA3P2', 'FGD3', 'RPL36AP30', 'RASGEF1C', 'BNIP2', 'ARF4', 'PRR4', 'TAS2R14', 'PRH1', 'FLT4', 'NELFB', 'AMPD3', 'DSP-AS1', 'RPS6KA1', 'IGKV1OR2-108', 'TTLL9', 'MRPL16', 'RIPK1', 'RFX1', 'PEBP1P1', 'PPIC-AS1', 'SMG1P3', 'SLC7A5P2', 'SMIM10L1', 'SNORA14A', 'TNS2-AS1', 'HNRNPA1P49', 'CASTOR3P', 'DMBX1', 'KCTD4', 'TM2D1', 'NCOR2', 'PRELID1', 'BAG3', 'ZNF232-AS1', 'NUP43', 'HK2', 'ABCD1P2', 'EMSLR', 'API5P2', 'RPS10P7', 'RN7SKP114', 'MFHAS1', 'AFDN-DT', 'ZYG11B', 'S100PBP', 'LINC01600', 'EIF3B', 'GPR161', 'H3-3A', 'LGALSL', 'PTDSS2', 'GNPTAB', 'PTPRVP', 'INTS6P1', 'PKN2-AS1', 'CSTB', 'RNU6-72P'), ],
             annotation_col = TCGA.annotations["Sample_Group"],
             cluster_distance_col = "correlation",
             cluster_distance_row = "correlation")




pheatmap::pheatmap(mat = data_Promoters[rownames(data_Promoters) %in%
                                        rownames(A1.now[A1.now$DM == TRUE,]), ],
                   annotation_col = TCGA.annotations[c("Sample_Group", "ER_status", "gender", "retrospective_collection")],
                   cluster_cols = T,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   scale = "column")


                   #cluster_distance_col = "correlation",
                   #cluster_distance_row = "correlation")





# Heatmap with differentially methylated genes with our samples and Australian ones

FC.threshold     = 1
p.val.threshold  = 0.05
data.pro.A1$"DifferentiallyMethylated" <- ifelse(data.pro.A1$adj.P.Val <= p.val.threshold & abs(data.pro.A1$logFC) >= FC.threshold, "TRUE", "FALSE")
data.pro.A2$"DifferentiallyMethylated" <- ifelse(data.pro.A2$adj.P.Val <= p.val.threshold & abs(data.pro.A2$logFC) >= FC.threshold, "TRUE", "FALSE")


heatmap_func(data = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% data.pro.A1[data.pro.A1$DifferentiallyMethylated == TRUE, ]$Gene, ],
             annotations       = annotations,
             subsetAnnotation  = rev(c("Study", "Subtype", "ER_status", "Variant_type", "ClinVar", "ATM_LOH")),
             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
             fontsize_row      = 8,
             fontsize_col      = 10,
             main              = 'Heatmap with promotors found differentially methylated when comparing tumours \n from pathogenic of likely pathogenic variant carriers to control tumours',
             )






Heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("HNRNPA1P49"),
                            annotations_A1$Sample_Name],
                        annotations       = annotations_A1,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")






# genes of random forest on our samples and Australian ones using ALL promoters

# A1
## With only A1 samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PDIA3P2", "CSTB", "SMG1P3", "RPS6KA1", "PRR4", "OR52H1", "MALT1-AS1", "POLR2L", "HNRNPA1P49"),
                            annotations_A1$Sample_Name],
                        annotations       = annotations_A1,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")

## With all samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PDIA3P2", "SCAPER", "PKN2-AS1", "INTS6P1", "TAS2R14", "POLR2L", "RASGEF1C", "CSTB"),
                            annotations$Sample_Name],
                        annotations       = annotations,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")






# A2
## With only A2 samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("CSTB", "PDIA3P2", "HNRNPA1P49", "OR52H1", "SMG1P3", "MALT1-AS1", "PRR4", "PRH1"),
                            annotations_A2$Sample_Name],
                        annotations       = annotations_A2,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")

## With all samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("CSTB", "PDIA3P2", "HNRNPA1P49", "OR52H1", "SMG1P3", "MALT1-AS1", "PRR4", "PRH1"),
                            annotations$Sample_Name],
                        annotations       = annotations,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")










# genes of random forest on our samples and Australian ones using DM promoters

# A1
## With only A1 samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PDIA3P2", "SCAPER", "PKN2-AS1", "INTS6P1", "TAS2R14", "POLR2L", "RASGEF1C", "CSTB"),
                            annotations_A1$Sample_Name],
                        annotations       = annotations_A1,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")

## With all samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PDIA3P2", "SCAPER", "PKN2-AS1", "INTS6P1", "TAS2R14", "POLR2L", "RASGEF1C", "CSTB"),
                            annotations$Sample_Name],
                        annotations       = annotations,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")


#A2
## With only A2 samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PKN2-AS1", "PDIA3P2", "SCAPER", "PRH1", "INTS6P1", "SLC7A5P2", "RASGEF1C", "MALT1-AS1"),
                            annotations_A2$Sample_Name],
                        annotations       = annotations_A1,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")

## With all samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PKN2-AS1", "PDIA3P2", "SCAPER", "PRH1", "INTS6P1", "SLC7A5P2", "RASGEF1C", "MALT1-AS1"),
                            annotations$Sample_Name],
                        annotations       = annotations,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")









heatmap <- heatmap_func(data              = data_Promoters[rownames(data_Promoters) %in% c("PKN2-AS1", "PDIA3P2", "SCAPER", "PRH1", "INTS6P1", "SLC7A5P2", "RASGEF1C", "MALT1-AS1"),
                            annotations$Sample_Name],
                        annotations       = TCGA.annotations,
                        #subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")


pheatmap::pheatmap(mat = data_Promoters[rownames(data_Promoters) %in% c("PKN2-AS1", "PDIA3P2", "SCAPER", "PRH1", "INTS6P1", "SLC7A5P2", "RASGEF1C", "MALT1-AS1"),],
                  annotation_col= TCGA.annotations[c("Sample_Group", "ER_status", "gender", "retrospective_collection")],
                                        #subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                  color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                  fontsize_row      = 8,
                  fontsize_col      = 10,
                  main              = "Heatmap")


























heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PDIA3P2", "CSTB", "SMG1P3", "RPS6KA1", "PRH1", "OR52H1", "MALT1-AS1", "POLR2L", "HNRNPA1P49"),
                            annotations_A1$Sample_Name],
                        annotations       = annotations_A1,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")



heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PDIA3P2", "CSTB", "SMG1P3", "PRH1", "OR52H1", "HNRNPA1P49", "PBRM1"),
                            annotations_A2$Sample_Name],
                        annotations       = annotations_A2,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type", "ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")




heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PDIA3P2", "CSTB", "SMG1P3", "RPS6KA1", "PRH1", "OR52H1", "MALT1-AS1", "POLR2L", "HNRNPA1P49"),
                            annotations$Sample_Name],
                        annotations       = annotations,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type", "ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")



heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PDIA3P2", "CSTB", "SMG1P3", "PRH1", "OR52H1", "HNRNPA1P49", "PBRM1"),
                            annotations$Sample_Name],
                        annotations       = annotations,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type", "ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")












ssgsea = function(X, genesets, alpha = 0.25, scale = T, norm = F, single = T) {
    rownames = rownames(X)
    numgenes = nrow(X)
    genesets = lapply(genesets, function(genes) {which(rownames %in% genes)})
    
            # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

            # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(Rcol) {
        generanks = order(Rcol, decreasing = TRUE)
        
                                # Calc es for each gene set
        essample = sapply(genesets, function(genesetidx) {
                                        # pos: match (within the gene set)
                                                # neg: non-match (outside the gene set)
            indicatorpos = generanks %in% genesetidx
            indicatorneg = !indicatorpos
            
            rankalpha  = (Rcol[generanks] * indicatorpos) ^ alpha
            
            stepcdfpos = cumsum(rankalpha)    / sum(rankalpha)
            stepcdfneg = cumsum(indicatorneg) / sum(indicatorneg)
            
            stepcdfdiff = stepcdfpos - stepcdfneg

                                        # Normalize by gene number
            if (scale) stepcdfdiff = stepcdfdiff / numgenes
            
                                        # Use ssGSEA or not
            if (single) {
                sum(stepcdfdiff)
            } else {
                stepcdfdiff[which.max(abs(stepcdfdiff))]
            }
        })
        unlist(essample)
    })
    
    if (length(genesets) == 1) es = matrix(es, nrow = 1)
    
                                        # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))
    
                                        # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}


data.frame(row.names = A1.Pro.KeggGse$Description)

colnames <- A1.Pro.KeggGse$Description
geneLists <- stringr::str_split(A1.Pro.KeggGse$core_enrichment, "/")

max = 0
for (i in 1:length(geneLists)){
     if ( length(geneLists[[i]]) > max){
         max = length(geneLists[[i]])
     }
}


res = data.frame(var=1:max)
for (i in 1:length(geneLists))
{
    print(geneLists[[i]])
    res[colnames[[i]]] <- c(geneLists[[i]], rep("", nrow(res)-length(geneLists[[i]])))
}


gene_sets = as.list(as.data.frame(res))


assign('resFinal', ssgsea(as.matrix(M.Promoters.annotated), gene_sets, scale = TRUE, norm = FALSE))
resFinal = resFinal[! rownames(resFinal) == "var",]
sauv = resFinal



#zscore the ssgsea output for comparative analysis
resFinal = (resFinal - rowMeans(resFinal))/(matrixStats::rowSds(as.matrix(resFinal)))[row(resFinal)]


identical(rownames(annotations), colnames(resFinal))


pheatmap::pheatmap(resFinal,
                   annotation_col = annotations[c("Sample_Group", "Variant_type_aggregated")])




resFinal = t(resFinal)





#mat = mat[, order(colnames(mat))]
#info = info[order(rownames(info)), ]








# New random forest
data = read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation/marqueur_identification_A1.csv")

gene = data[data$Count >= 50, ]$Gene

# genes of random forest on our samples and Australian ones using ALL promoters

# A1
## With only A1 samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% gene,
                            annotations_A1$Sample_Name],
                        annotations       = annotations_A1,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")

## With all samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("PDIA3P2", "SCAPER", "PKN2-AS1", "INTS6P1", "TAS2R14", "POLR2L", "RASGEF1C", "CSTB"),
                            annotations$Sample_Name],
                        annotations       = annotations,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")






# A2
## With only A2 samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("CSTB", "PDIA3P2", "HNRNPA1P49", "OR52H1", "SMG1P3", "MALT1-AS1", "PRR4", "PRH1"),
                            annotations_A2$Sample_Name],
                        annotations       = annotations_A2,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")

## With all samples
heatmap <- heatmap_func(data              = M.Promoters.annotated[rownames(M.Promoters.annotated) %in% c("CSTB", "PDIA3P2", "HNRNPA1P49", "OR52H1", "SMG1P3", "MALT1-AS1", "PRR4", "PRH1"),
                            annotations$Sample_Name],
                        annotations       = annotations,
                        subsetAnnotation  = c("ATM_LOH", "ClinVar", "Variant_type","ER_status",  "Study"),
                        color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                        fontsize_row      = 8,
                        fontsize_col      = 10,
                        main              = "Heatmap")











p.val.threshold = 0.05
FC.threshold = 1

d1 <- data.pro.A1
d1$DM <- ifelse(d1$adj.P.Val <= p.val.threshold & abs(d1$logFC) >= FC.threshold, "TRUE", "FALSE")
d1$comp <- "PV"

d2 <- data.pro.AVUS
d2$DM <- ifelse(d2$adj.P.Val <= p.val.threshold & abs(d2$logFC) >= FC.threshold, "TRUE", "FALSE")
d2$comp <- "VUS"

d3 <- data.pro.A0
d3$DM <- ifelse(d3$adj.P.Val <= p.val.threshold & abs(d3$logFC) >= FC.threshold, "TRUE", "FALSE")
d3$comp <- "All"

p <- rbind(d1, d2, d3)

#p <- merge(d1, d2, by = "Gene")
#p <- merge(, data.pro.AVUS, by = "Gene")


p <- p[order(p$Gene),]

DM.PV <- p[p$comp == "PV" & p$DM == "TRUE",]$Gene
DM.VUS <- p[p$comp == "VUS" & p$DM == "TRUE",]$Gene
DM.All <- p[p$comp == "All" & p$DM == "TRUE",]$Gene


# Are data normally distributed ?
shapiro.test(p[p$comp == "PV",]$logFC)

ks.test(x=p[p$comp == "PV",]$logFC,y='pnorm',alternative='two.sided')






# comparison PV vs All and VUS vs All
p$col.PV.VUS <- ifelse(p$Gene %in% DM.PV & p$Gene %in% DM.VUS, "orange",
                       ifelse(p$Gene %in% DM.PV & ! p$Gene %in% DM.VUS, "red",
                              ifelse(p$Gene %in% DM.VUS & ! p$Gene %in% DM.PV, "yellow", "black"))
                       )

a <- ggplot() + geom_point(aes(x=p$logFC[p$comp == "PV"],
                               y=p$logFC[p$comp == "VUS"]), 
                           color = p[p$comp %in% c("VUS", "PV"),c("Gene", "col.PV.VUS")][!duplicated(p[p$comp %in% c("VUS", "PV"),c("Gene", "col.PV.VUS")]),]$col.PV.VUS) +
    xlab("logFC of ATM PV vs All controls") + ylab("logFC of ATM VUS vs All controls") + 
    theme_bw()




# comparison PV vs All and All vs All
p$col.PV.All <- ifelse(p$Gene %in% DM.PV & p$Gene %in% DM.All, "orange",
                       ifelse(p$Gene %in% DM.PV & ! p$Gene %in% DM.All, "red",
                              ifelse(p$Gene %in% DM.All & ! p$Gene %in% DM.PV, "yellow", "black"))
                       )

b <- ggplot() + geom_point(aes(x=p$logFC[p$comp == "PV"],
                          y=p$logFC[p$comp == "All"]), 
                          color = p[p$comp %in% c("All", "PV"),c("Gene", "col.PV.All")][!duplicated(p[p$comp %in% c("All", "PV"),c("Gene", "col.PV.All")]),]$col.PV.All) +
    xlab("logFC of ATM PV vs All controls") + ylab("logFC of all ATM vs all controls") + 
    theme_bw()






# comparison VUS vs All and All vs All
p$col.VUS.All <- ifelse(p$Gene %in% DM.VUS & p$Gene %in% DM.All, "orange",
                       ifelse(p$Gene %in% DM.VUS & ! p$Gene %in% DM.All, "red",
                              ifelse(p$Gene %in% DM.All & ! p$Gene %in% DM.VUS, "yellow", "black"))
                       )

c <- ggplot() + geom_point(aes(x=p$logFC[p$comp == "VUS"],
                          y=p$logFC[p$comp == "All"]), 
                          color = p[p$comp %in% c("All", "VUS"),c("Gene", "col.VUS.All")][!duplicated(p[p$comp %in% c("All", "VUS"),c("Gene", "col.VUS.All")]),]$col.VUS.All) +
    xlab("logFC of ATM VUS vs All controls") + ylab("logFC of all ATM vs All controls") + 
    theme_bw()




ggpubr::ggarrange(a, b, c,
                  labels = c("A", "B", "C"))









# Load data
OutFiles = file.path(OutFiles1, "Kegg/Gse/")
load(file.path(OutFiles, "AllKeggGseObjects.Rdata"))

load(file.path(Origin, "/matrices/Mapping.RData"))


A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "Oxidative phosphorylation", ]


HeatmapNew <- function(data,
                       ListGenes,
                       annotations = annotations_A0A2,
                       show_colnames = TRUE,
                       scale = "none",
                       main = "")
{
    newnames <- lapply(
        rownames(data[ListGenes, rownames(annotations)]),
        function(x) bquote(italic(.(x))))
    #print(newnames)
    heat_annotations <- data.frame(row.names     = rownames(annotations),
                                   Study         = annotations$Study,
                                   Sample_Group  = annotations[['Sample_Group']],
                                   ER_status     = as.character(annotations[['ER_status']]),
                                   Variant_type  = annotations$Variant_type_aggregated,
                                   ClinVar       = annotations$ClinVar
                                   )
    print(head(heat_annotations))
    heat_annotations$ClinVar <- ifelse(heat_annotations$ClinVar %in% c("Pathogenic/Likely pathogenic", "Pathogenic", "Likely pathogenic"), "PV", heat_annotations$ClinVar)
    heat_annotations$ClinVar <- ifelse(heat_annotations$ClinVar %in% c("Uncertain significance"), "VUS", heat_annotations$ClinVar)
    heat_annotations$ClinVar <- ifelse(heat_annotations$ClinVar == "", NA, heat_annotations$ClinVar)
    heat_annotations$Variant_type <- ifelse(heat_annotations$Variant_type == "", NA, heat_annotations$Variant_type)
    heat_annotations$Sample_Group <- ifelse(heat_annotations$Sample_Group == "Non.ATM", "Control", heat_annotations$Sample_Group) 
    anno.col = list(
        ClinVar = c(PV = "#E199FF",
                    VUS = "#FBA625"),
        Variant_type = c(LoF = "#95CA00",
                         MV = "#01DAE1"),
        ER_status = c("0" = "#FF82CB",
                      "1" = "#00DCA6"),
        Sample_Group = c(ATM = "#00DB97",
                         Control = "#B0ABFF"),
        Study = c(ABCFR = "#FF9389",
                  CoFAT = "#82B6FF",
                  GENESIS = "#D3B900",
                  MCCS = "#FF83FF")
    )
    data <- data[ListGenes, rownames(annotations)]
  #  data <- data[rowSums(is.na(data)) != ncol(data), ]
    plot <- pheatmap(data,
                     annotation_col = heat_annotations,
                     annotation_colors = anno.col,
                     color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                     scale = scale,
                     clustering_distance_cols = "correlation",
                     labels_row = as.expression(newnames),
                     show_colnames = show_colnames,
                     main = main)
    
    return(plot)
}


List <- ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "Hedgehog signaling pathway", ],
                           column = "core_enrichment")
List <- List$core_enrichment

HeatmapNew(data = data_Promoters,
           ListGenes = List,
           annotations = annotations_A2)


# Try to predict with those biomarkers
data <- merge(t(data_Promoters[, rownames(annotations_A2)]), annotations_A2["Sample_Group"], by = "row.names", all.x=T, all.y=F)
data <- data[c(List, "Sample_Group")]

names(data) <- make.names(names(data))
data$Sample_Group <- factor(data$Sample_Group)


Acc <- c()
F1 <- c()
for (i in c(1:500))
{
# Split train/test
    sample <- caTools::sample.split(data$Sample_Group, SplitRatio = 0.7)
    train  <- subset(data, sample == TRUE)
    test   <- subset(data, sample == FALSE)
    print(table(train$Sample_Group))
    print(table(test$Sample_Group))
# Random forest
    rf1 <- randomForest::randomForest(Sample_Group ~ ., data = train,
                                      ntree = 1000, mtry = dim(data)[2],
                                      importance=TRUE)
# Predict
    pred <- predict(rf1, newdata = test, type = "response")
    cm <- table(test$Sample_Group,pred)
    cm
# Confusion matrix
    res1 <- caret::confusionMatrix(cm, mode = "everything")
#    print(res1$overall[1])
# Add accuracy
    Acc <- c(Acc, res1$overall[["Accuracy"]])
    F1 <- c(F1, res1$byClass[["F1"]])
}

min(Acc)
max(Acc)
mean(Acc)
median(Acc)


min(F1)
max(F1)
mean(F1)
median(F1)

plot(hist(F1))

summary(F1)



#the train AUC
rf_p_train <- predict(rf1, type="prob")[,2]
rf_pr_train <- ROCR::prediction(rf_p_train, train$Sample_Group)
r_auc_train <- ROCR::performance(rf_pr_train, measure = "auc")@y.values[[1]]
r_auc_train   #0.8098

#the test AUC
rf_p_test <- predict(rf1, type="prob",newdata = test)[,2]
rf_pr_test <- ROCR::prediction(rf_p_test, test$Sample_Group)
r_auc_test <- ROCR::performance(rf_pr_test, measure = "auc")@y.values[[1]]
r_auc_test    #0.956



# Other method: get all genes from core enrichment, get their MeanDecreaseAccuracy and make a clustering with them
ScoreGenes <- function(matrix,
                         EnrichmentResult,
                         anno)
{
    res = list()
    for (n in EnrichmentResult$Description)
    {
        #print(n)     
        List <- ParseEnrichResults(enrichResultObject = EnrichmentResult[EnrichmentResult$Description == n, ],
                                   column = "core_enrichment")
        List <- unlist(List$core_enrichment)
        print(List)
        res <- append(res, List)
    }
    res <- unlist(res)
    res <- List[!duplicated(res)]
    List <- res
    
    # Try to predict with those biomarkers  
    data <- merge(t(matrix[, rownames(anno)]), anno["Sample_Group"], by = "row.names", all.x=T, all.y=F)
#    print(head(data))
    data <- data[, c(List, "Sample_Group")]
    
    ## Set seed for reproducibility
    set.seed(123)

    print("step 1")
    ## Define repeated cross validation with 5 folds and three repeats
    repeat_cv <- caret::trainControl(method='repeatedcv', number=5, repeats=3)
    print("step 2")
    ## Split the data so that we use 70% of it for training
    train_index <- createDataPartition(y=data$Sample_Group, p=0.7, list=FALSE)
    print("step 3")
    ## Subset the data
    training_set <- iris[train_index, ]
    testing_set <- iris[-train_index, ]
    
    return(res)
}

Data1 <- ScoreGenes(matrix = data_Promoters,
                     EnrichmentResult = A2.Pro.KeggGse,
                     anno = annotations_A2)





# For each enriched pathways test the accuracy of the classification
testPathways <- function(matrix,
                         EnrichmentResult,
                         anno)
{
    res = list()
    for (n in EnrichmentResult$Description)
    {
        print(n)     
        List <- ParseEnrichResults(enrichResultObject = EnrichmentResult[EnrichmentResult$Description == n, ],
                                   column = "core_enrichment")
        List <- List$core_enrichment
        
        # Try to predict with those biomarkers  
        data <- merge(t(matrix[, rownames(anno)]), anno["Sample_Group"], by = "row.names", all.x=T, all.y=F)
#        rownames(data) <- data$Row.names
        data <- data[c(List, "Sample_Group")]
#        print(head(data))
        names(data) <- make.names(names(data))
        data$Sample_Group <- factor(data$Sample_Group)
        # empty vetcor
        Acc <- c()
        F1 <- c()
        #print(head(data))
        for (i in c(1:10))
        {
            # Split train/test
            sample <- caTools::sample.split(data$Sample_Group, SplitRatio = 0.7)
            train  <- subset(data, sample == TRUE)
            test   <- subset(data, sample == FALSE)
            #print(table(train$Sample_Group))
            #print(table(test$Sample_Group))
            # Random forest
            rf1 <- randomForest::randomForest(Sample_Group ~ ., data = train,
                                              ntree = 1000, mtry = dim(data)[2]/3,
                                              importance=TRUE)
            imp <- randomForest::importance(rf1, type=1, scale = F) # permutation importances
            print(imp)
            #print(rf1$importance[order(rf1$importance[, 1], decreasing = TRUE), ])
            # Predict
            pred <- predict(rf1, newdata = test, type = "response")
            cm <- table(test$Sample_Group,pred)
            cm
            # Confusion matrix
            res1 <- caret::confusionMatrix(cm, mode = "everything")
            # Add accuracy
            Acc <- c(Acc, round(res1$overall[["Accuracy"]], 3))
            F1 <- c(F1, round(res1$byClass[["F1"]], 3))
        }
        #the train AUC
        rf_p_train <- predict(rf1, type="prob")[,2]
        rf_pr_train <- ROCR::prediction(rf_p_train, factor(train$Sample_Group))
        r_auc_train <- ROCR::performance(rf_pr_train, measure = "auc")@y.values[[1]]
        r_auc_train   #0.8098
        #the test AUC
        rf_p_test <- predict(rf1, type="prob",newdata = test)[,2]
        rf_pr_test <- ROCR::prediction(rf_p_test, factor(test$Sample_Group))
        r_auc_test <- ROCR::performance(rf_pr_test, measure = "auc")@y.values[[1]]
        r_auc_test    #0.956
        res[[n]] <- list("name" = n,
                         "Accuracy" = Acc,
                         "F1" = F1,
                         "AUC_train" = r_auc_train,
                         "AUC_test" = r_auc_test)
    }

    columns = c("Pathway", "Accuracy", "Summary_Acc",
                "F1",
                "Summary_F1",
                "AUC_train", "AUC_test")
    Data = data.frame(matrix(nrow = length(res), ncol = length(columns)))
    colnames(Data) = columns
    
    for (n in 1:length(res))
    {
        Data[n,] <- c(res[[n]]$name ,
                      paste(list(res[[n]]$Accuracy), collapse=', '), 
                      paste(list(summary(res[[n]]$Accuracy)), collapse=', '), 
                      paste(list(res[[n]]$F1), collapse=', '), 
                      paste(list(summary(res[[n]]$F1)), collapse=", "), 
                      res[[n]]$AUC_train,
                      res[[n]]$AUC_test)
    }
    
    return(Data)
}

Data <- testPathways(matrix = data_Promoters,
                     EnrichmentResult = A2.Pro.KeggGse,
                     anno = annotations_A2)

write.csv(Data, "~/ATM_Analysis/svn/analyse/results/FinalMethylation/test_biomarkers_A2_v2.2.csv")


Data <- testPathways(matrix = data_Promoters,
                     EnrichmentResult = A0A2.Pro.KeggGse,
                     anno = annotations_A0A2)

write.csv(Data, "~/ATM_Analysis/svn/analyse/results/FinalMethylation/test_biomarkers_A0A2.2.csv")


Data <- testPathways(matrix = data_Promoters,
                     EnrichmentResult = AVUSA2.Pro.KeggGse,
                     anno = annotations_AVUSA2)

write.csv(Data, "~/ATM_Analysis/svn/analyse/results/FinalMethylation/test_biomarkers_AVUSA2.2.csv")








List = c("GPR161", "MEGF8", "GSK3B", "SMO", "EVC", 
"CSNK1G2", "SMURF2", "MOSMO", "EFCAB7", "SUFU", "GAS1", "CUL1", "LRP2", "EVC2", "ARRB1", "GLI3", "CSNK1G1"
)










List <- ParseEnrichResults(enrichResultObject = A0A2.Pro.KeggGse[A0A2.Pro.KeggGse$Description == "Pathways in cancer", ],
                           column = "core_enrichment")
List <- List$core_enrichment

List <- ParseEnrichResults(enrichResultObject = A0A2.Pro.KeggGse[A0A2.Pro.KeggGse$Description == "Taste transduction", ],
                           column = "core_enrichment")
List <- List$core_enrichment


HeatmapNew(data = data_Promoters,
           ListGenes = List,
           annotations = annotations_A0A2)


pathways <- c("OXIDATIVE PHOSPHORYLATION",
"DIABETIC CARDIOMYOPATHY", 
"PROTEASOME", 
"CHEMICAL CARCINOGENESIS - REACTIVE OXYGEN SPECIES", 
"PRION DISEASE", 
"HUNTINGTON DISEASE", 
"AMYOTROPHIC LATERAL SCLEROSIS", 
"ALZHEIMER DISEASE", 
"THERMOGENESIS", 
"NON-ALCOHOLIC FATTY LIVER DISEASE", 
"PARKINSON DISEASE", 
"SPINOCEREBELLAR ATAXIA", 
"PATHWAYS OF NEURODEGENERATION - MULTIPLE DISEASES")


pathways <- c"CHEMICAL CARCINOGENESIS - DNA ADDUCTS", 
"CAFFEINE METABOLISM", 
"DRUG METABOLISM - CYTOCHROME P450", 
"STEROID HORMONE BIOSYNTHESIS", 
"RETINOL METABOLISM", 
"METABOLISM OF XENOBIOTICS BY CYTOCHROME P450")


pathways <- c("GLIOMA", 
"BREAST CANCER", 
"HUMAN PAPILLOMAVIRUS INFECTION", 
"HUMAN T-CELL LEUKEMIA VIRUS 1 INFECTION", 
"HEPATOCELLULAR CARCINOMA", 
"PATHWAYS IN CANCER", 
"ENDOCRINE RESISTANCE", 
"GASTRIC CANCER", 
"CHRONIC MYELOID LEUKEMIA")

pathways <- str_to_sentence(pathways)

List <- ParseEnrichResults(enrichResultObject = A0A2.Pro.KeggGse[A0A2.Pro.KeggGse$Description %in% pathways, ],
                           column = "core_enrichment")
List <- List$core_enrichment
List <- List[!duplicated(List)]

HeatmapNew(data = data_Promoters_Beta,
           ListGenes = List,
           annotations = annotations_A0A2,
           scale = "none")


data <- merge(t(data_Promoters[, rownames(annotations_A0A2)]), annotations_A0A2["Sample_Group"], by = "row.names", all.x=T, all.y=F)
#data <- data[c(List, "Sample_Group")]

names(data) <- make.names(names(data))
data$Sample_Group <- factor(data$Sample_Group)

# options(expressions = 5e5)


res <- list()
for (n in c(1:10))
    {
        fit <- randomForest::randomForest(Sample_Group ~ ., data = data,
                                          ntree = 300)
        if (n==1)
        {
            res <- data.frame(fit$importance)
        } else{
            res <- merge(res, data.frame(fit$importance), by="row.names")
            rownames(res) <- res$Row.names
            res <- subset(res, select= - Row.names)
        }
    }



fit <- randomForest::randomForest(Sample_Group ~ ., data = data,
                                  ntree = 5000)


randomForest::varImpPlot(fit)

eeeeeeeeeeeeeeee



#pheatmap(data_Promoters[List, rownames(annotations_A0A2)],
#         annotation_col = heat_annotations,
#         scale = 'column')

#pheatmap(data_Promoters_Beta[List, rownames(annotations_A0A2)],
#         annotation_col = heat_annotations,
#         scale = 'none')





pheatmap(data_Promoters_Beta[pro.A0A2$Gene,rownames(annotations_A0A2)],
         annotation_col = heat_annotations,
         scale = 'none')



    
plotByGene(M.Promoters.annotated,
           Gene = "ATM",
           Annotations = annotations_A0A2)



plotByGene(M.Promoters.annotated,
           Gene = "ATM",
           Annotations = annotations_A0A2,
           facet = "ClinVar") 

annotations_A0A2[["ClinVar"]]



plotByGene(M.Promoters.annotated,
           Gene = "GPR161",
           Annotations = annotations_A0A2)
    

plotByGene(M.Promoters.annotated,
           Gene = "CDK13",
           Annotations = annotations_A0A2)




plot <- function(data = M.Promoters.annotated,
                 Annotations = annotations_A0A2,
                 Gene = "ATM")
{
    require(ggpubr)
    require(dplyr)

    Annotations$ATM_LOH  <- ifelse(Annotations$Sample_Name %in% ATM.LOH.Yes, "LOH",
                            ifelse(Annotations$Sample_Name %in% ATM.LOH.No, "No LOH",
                                   "Unknown LOH")
                            )
    a <- data.frame(row.names    = rownames(Annotations),
                    Sample_Group = Annotations[["Sample_Group"]],
                    ClinVar      = Annotations[["ClinVar"]],
                    ATM_LOH      = Annotations[["ATM_LOH"]])
    a[a$ClinVar %in% c( "Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"),]$ClinVar <- "ATM PV ER+"
    print(head(a))
    print(a$ClinVar)
    if ("Uncertain significance" %in% a$ClinVar)
    {
        a[a$ClinVar %in% c("Uncertain significance"),]$ClinVar <- "ATM VUS ER+"
    }
    if ("" %in% a$ClinVar)
    {
        a[a$ClinVar %in% c(""),]$ClinVar <- "Sporadic"
    }
    print(a[a$ClinVar %in% c(""),]$ClinVar)
    a$ClinVar <- factor(a$ClinVar, levels = c("ATM PV ER+", "ATM VUS ER+", "Sporadic"))
    a <- merge(a, t(data[Gene, ]), by = "row.names")
    names(a) <- c("Sample", "Sample_Group", "ClinVar", "ATM_LOH", "Beta")

    a$ATM_LOH <- ifelse(a$Sample_Group == "Non.ATM", "Sporadic",
                        paste0(a$Sample_Group, " ", a$ATM_LOH))
    a$ATM_LOH <- factor(a$ATM_LOH, levels = c("ATM LOH", "ATM No LOH", "ATM Unknown LOH", "Sporadic"))
    print(a)
    p1 <- ggviolin(a,
               x        = "ATM_LOH",
               y        = "Beta",
               palette  = c("#B72B22", "#FF5349", "#F8756D", "#01BEC5"),
               fill     = "ATM_LOH",
               add      = "boxplot") +
    theme(legend.position="none") +#,
#          plot.title = element_text(hjust = 0.5),
#          plot.subtitle = element_text(hjust = 0.5, size=rel(0.95))) + 
                                        #stat_compare_means(comparisons  = comparisons, method="t.test")# +
#        stat_pvalue_manual(t.test, label = "label", y.position = max) +
        ylab("M values") +
            xlab("Type of variant with LOH")
#        custom.theme
    return(p1)
}



plot(data = M.Promoters.annotated,
     Annotations = annotations_A0A2,
     Gene = "ATM")


p1 <- ggpubr::ggviolin(M.Promoters.annotated["ATM", ],
               x        = "Sample_Group",
               y        = "mean",
               palette  = col, #c("#01BEC5", "#F8756E"),
               fill     = "Sample_Group",
               add      = "boxplot") +
    labs(title = paste0(main),
         subtitle = paste0("(", nprobes, " probes)")) +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size=rel(0.95))) + 
                                        #stat_compare_means(comparisons  = comparisons, method="t.test")# +
        stat_pvalue_manual(t.test, label = "label", y.position = max) +
        xlab("Type of tumours") +
        custom.theme






# Random Forest from Python

#A0A2
p <- HeatmapNew(data = data_Promoters,
           ListGenes = c("RFX1", "UBE2R2", "AMPD3", "INTS6P1", "AIMP2", "DNAJC7", "PABIR1", "SNORD88C", "SMCO3"),
           annotations = annotations_A0A2,
           show_colnames = FALSE,
           scale = "column")

ggsave(file.path(OutFiles, "A0A2signature_promoters.png"), p, width=400, height=275, units="mm")

p <- HeatmapNew(data = data_Promoters,
           ListGenes = c("RFX1", "UBE2R2", "AMPD3", "INTS6P1", "AIMP2", "DNAJC7", "PABIR1", "SNORD88C", "SMCO3"),
           annotations = annotations,
           show_colnames = FALSE,
           scale = "column")

ggsave(file.path(OutFiles, "A0A2signature_promoters_allSamplesWithERminus.png"), p, width=400, height=275, units="mm")


#A2
p <- HeatmapNew(data = data_Promoters,
           ListGenes = c("RFX1", "AMPD3", "INTS6P1", "HSPE1P10", "SNORD88C", "RIPK1", "MALT1-AS1", "AIMP2", "SUMO2P15", "UBE2R2", "PABIR1"),
           annotations = annotations_A2,
           show_colnames = FALSE,
           scale="column")

ggsave(file.path(OutFiles, "A2signature_promoters.png"), p, width=400, height=275, units="mm")

#signature A2 appliquée à tous les échantillons ER+ => fonctionne très bien
p <- HeatmapNew(data = data_Promoters,
           ListGenes = c("RFX1", "AMPD3", "INTS6P1", "HSPE1P10", "SNORD88C", "RIPK1", "MALT1-AS1", "AIMP2", "SUMO2P15", "UBE2R2", "PABIR1"),
           annotations = annotations_A0A2,
           show_colnames = FALSE,
           scale="column")

ggsave(file.path(OutFiles, "A2signature_promoters_allSamples.png"), p, width=400, height=275, units="mm")

#signature A2 appliquée à tous les échantillons ER+ et ER- 
p <- HeatmapNew(data = data_Promoters,
           ListGenes = c("RFX1", "AMPD3", "INTS6P1", "HSPE1P10", "SNORD88C", "RIPK1", "MALT1-AS1", "AIMP2", "SUMO2P15", "UBE2R2", "PABIR1"),
           annotations = annotations,
           show_colnames = FALSE,
           scale="column")

ggsave(file.path(OutFiles, "A2signature_promoters_allSamplesWithERminus.png"), p, width=400, height=275, units="mm")














data_Promoters.TCGA <- read.csv("~/ATM_Analysis/svn/analyse/results/FinalMethylation_TCGA/matrices/M_Promoters_annotated.csv", row.names=1)

#A0A2
HeatmapNew(data = data_Promoters.TCGA,
           ListGenes = c("RFX1", "UBE2R2", "AMPD3", "INTS6P1", "AIMP2", "DNAJC7", "PABIR1", "SNORD88C", "SMCO3"),
           annotations = annotations_A0A2,
           show_colnames = FALSE,
           scale = "column")


#A2
HeatmapNew(data = data_Promoters.TCGA,
           ListGenes = c("RFX1", "AMPD3", "INTS6P1", "HSPE1P10", "SNORD88C", "RIPK1", "MALT1-AS1", "AIMP2", "SUMO2P15", "UBE2R2", "PABIR1"),
           annotations = annotations_A2,
           show_colnames = FALSE,
           scale="column")

#signature A2 appliquée à tous les échantillons ER+ => fonctionne très bien
HeatmapNew(data = data_Promoters.TCGA,
           ListGenes = c("RFX1", "AMPD3", "INTS6P1", "HSPE1P10", "SNORD88C", "RIPK1", "MALT1-AS1", "AIMP2", "SUMO2P15", "UBE2R2", "PABIR1"),
           annotations = annotations_A0A2,
           show_colnames = FALSE,
           scale="column")


data_Promoters.TCGA[c("RFX1", "AMPD3", "INTS6P1", "HSPE1P10", "SNORD88C", "RIPK1", "MALT1-AS1", "AIMP2", "SUMO2P15", "UBE2R2", "PABIR1"), rownames(annotations_A0A2)]







# From Pathways

# A2
p <- HeatmapNew(data = data_Promoters,
           ListGenes = ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "Taste transduction", ],
                                          column = "core_enrichment")$core_enrichment,
           annotations = annotations_A2,
           show_colnames = FALSE,
           scale="column",
           main = "Taste transduction")

ggsave(file.path(OutFiles, "A2_Taste_transduction_promoters.png"), p, width=400, height=275, units="mm")


p <- HeatmapNew(data = data_Promoters,
           ListGenes = ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "Hedgehog signaling pathway", ],
                                          column = "core_enrichment")$core_enrichment,
           annotations = annotations_A2,
           show_colnames = FALSE,
           scale="column",
           main = "Hedgehog signaling pathway")

ggsave(file.path(OutFiles, "A2_Hedgehog_signaling_pathway_promoters.png"), p, width=400, height=275, units="mm")


p <- HeatmapNew(data = data_Promoters,
           ListGenes = ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "Hedgehog signaling pathway", ],
                                          column = "core_enrichment")$core_enrichment,
           annotations = annotations,
           show_colnames = FALSE,
           scale="column",
           main = "Hedgehog signaling pathway")

ggsave(file.path(OutFiles, "A2_Hedgehog_signaling_pathway_promoters_allSamplesWithERminus.png"), p, width=400, height=275, units="mm")



p <- HeatmapNew(data = data_Promoters,
           ListGenes = ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "Ubiquitin mediated proteolysis", ],
                                          column = "core_enrichment")$core_enrichment,
           annotations = annotations_A2,
           show_colnames = FALSE,
           scale="column",
           main = "Ubiquitin mediated proteolysis")

ggsave(file.path(OutFiles, "A2_Ubiquitin_mediated_proteolysis_promoters.png"), p, width=400, height=275, units="mm")


p <- HeatmapNew(data = data_Promoters,
           ListGenes = ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "Ubiquitin mediated proteolysis", ],
                                          column = "core_enrichment")$core_enrichment,
           annotations = annotations,
           show_colnames = FALSE,
           scale="column",
           main = "Ubiquitin mediated proteolysis")

ggsave(file.path(OutFiles, "A2_Ubiquitin_mediated_proteolysis_promoters_allSamplesWithERminus.png"), p, width=400, height=275, units="mm")





p <- HeatmapNew(data = data_Promoters,
           ListGenes = ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "Alcoholism", ],
                                          column = "core_enrichment")$core_enrichment,
           annotations = annotations_A2,
           show_colnames = FALSE,
           scale="column",
           main = "Alcoholism")

ggsave(file.path(OutFiles, "A2_Alcoholism_promoters.png"), p, width=400, height=275, units="mm")



#p <- HeatmapNew(data = data_Promoters,
#           ListGenes = ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "Alcoholism", ],
#                                          column = "core_enrichment")$core_enrichment,
#           annotations = annotations_A2,
#           scale="column")
#ggsave(file.path(OutFiles, "A2_Alcoholism_promoters.png"), p)

#HeatmapNew(data = data_Promoters,
#           ListGenes = ParseEnrichResults(enrichResultObject = A2.Pro.KeggGse[A2.Pro.KeggGse$Description == "ATP-dependent chromatin remodeling", ],
#                                          column = "core_enrichment")$core_enrichment,
#           annotations = annotations_A2,
#           scale="column")



# A0A2
listPathways <- c("Olfactory transduction", "ATP-dependent chromatin remodeling", "Taste transduction", "Amyotrophic lateral sclerosis", "Hedgehog signaling pathway", "Hepatocellular carcinoma", "Ubiquitin mediated proteolysis")

for (i in listPathways)
{
    p <- HeatmapNew(data = data_Promoters,
                    ListGenes = ParseEnrichResults(enrichResultObject = A0A2.Pro.KeggGse[A0A2.Pro.KeggGse$Description == i, ],
                                                   column = "core_enrichment")$core_enrichment,
                    annotations = annotations_A0A2,
                    show_colnames = FALSE,
                    scale="column",
                    main = i)
    ggsave(file.path(OutFiles, paste0("A0A2_", i, "_promoters.png")), p, width=400, height=275, units="mm")
}


p <- HeatmapNew(data = data_Promoters,
           ListGenes = ParseEnrichResults(enrichResultObject = A0A2.Pro.KeggGse[A0A2.Pro.KeggGse$Description == "Hedgehog signaling pathway", ],
                                          column = "core_enrichment")$core_enrichment,
           annotations = annotations,
           show_colnames = FALSE,
           scale="column",
           main = "Hedgehog signaling pathway")

ggsave(file.path(OutFiles, "A0A2_Hedgehog_signaling_pathway_promoters_allSamplesWithERminus.png"), p, width=400, height=275, units="mm")



# A0AVUS
listPathways <- c("Cell cycle", 
                  "Autophagy - animal", 
                  "Spinocerebellar ataxia", 
                  "Salmonella infection", 
                  "Ribosome", 
                  "Pyrimidine metabolism", 
                  "Neurotrophin signaling pathway", 
                  "Central carbon metabolism in cancer", 
                  "Complement and coagulation cascades", 
                  "Signaling pathways regulating pluripotency of stem cells", 
                  "Gastric cancer", 
                  "Apoptosis", 
                  "Coronavirus disease - COVID-19", 
                  "Lipid and atherosclerosis", 
                  "Thermogenesis", 
                  "Natural killer cell mediated cytotoxicity")


for (i in listPathways)
{
    p <- HeatmapNew(data = data_Promoters,
                    ListGenes = ParseEnrichResults(enrichResultObject = AVUSA2.Pro.KeggGse[AVUSA2.Pro.KeggGse$Description == i, ],
                                                   column = "core_enrichment")$core_enrichment,
                    annotations = annotations_A0A2,
                    show_colnames = FALSE,
                    scale="column",
                    main = i)
    ggsave(file.path(OutFiles, paste0("AVUSA2_", i, "_promoters.png")), p, width=400, height=275, units="mm")
}





p.val.threshold = 0.05
FC.threshold = 1

d1 <- data.pro.A2
d1$DM <- ifelse(d1$adj.P.Val <= p.val.threshold & abs(d1$logFC) >= FC.threshold, "TRUE", "FALSE")
d1$comp <- "PV"

HeatmapNew(data = data_Promoters,
           ListGenes = d1[d1$DM == TRUE, ]$Gene,
           annotations = annotations_A2,
           show_colnames = FALSE,
           scale="column",
           main = "Heatmap of DM promoters with all ATM tumours")

















stop()

library(ggplot2)  # Load ggplot2 for violin plots
library(readxl)
library(UpSetR)

path <- file.path(Origin, "Classifier/final_analyses/A2_other_tests_validation_1Final_features_selected.xlsx")
exp1 <- read_excel(path)
colnames(exp1) <- paste0("exp1_", colnames(exp1), sep="") 

path <- file.path(Origin, "Classifier/final_analyses/A2_other_tests_validation_2Final_features_selected.xlsx")
exp2 <- read_excel(path)
colnames(exp2) <- paste0("exp2_", colnames(exp2), sep="") 


path <- file.path(Origin, "Classifier/final_analyses/A2_other_tests_validation_3Final_features_selected.xlsx")
exp3 <- read_excel(path)
colnames(exp3) <- paste0("exp3_", colnames(exp3), sep="") 


all.features <- c(as.list(exp1), as.list(exp2), as.list(exp3))
all.features <- lapply(all.features, function(x) x[!is.na(x)])






upset <- ggplotify::as.ggplot(upset(fromList(all.features),# sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"),
#      sets = rev(c("A2_Logistic Regression", "A2Inactive_Logistic Regression", "A0A2_Logistic Regression",
#                   "A2_Random Forest", "A2Inactive_Random Forest", "A0A2_Random Forest",
#                   "A2_XGBoost", "A2Inactive_XGBoost", "A0A2_XGBoost")),
                                    sets = rev(names(all.features)),
      keep.order = T, 
      empty.intersections = NULL,
      nsets = 9,
      sets.bar.color = "#56B4E9",
      order.by = "freq",
      text.scale = 2)
)

ggsave(filename = file.path(Origin, "/Classifier/upset.png"), upset,
       width = 50,
       height = 30,
       units = "cm",
       dpi=300)
                         


# Here I make the heatmaps


genes <- unique(c(all.features$"exp1_Logistic Regression",
                  all.features$"exp2_Logistic Regression",
                  all.features$"exp3_Logistic Regression"))


LR.1 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A2$Sample_Name],
                   annotations = annotations_A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames=T,
                   colAnno = c("Supposed biallelic inactivation", "Sample_Group"),
                   )

LR.2 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames=T,
                   colAnno = c("Supposed biallelic inactivation", "Sample_Group"),
                   )


genes <- unique(c(all.features$"exp1_Random Forest",
                  all.features$"exp2_Random Forest",
                  all.features$"exp3_Random Forest"))


RF.1 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A2$Sample_Name],
                   annotations = annotations_A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames=T,
                   colAnno = c("Supposed biallelic inactivation", "Sample_Group"),
                   )

RF.2 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames=T,
                   colAnno = c("Supposed biallelic inactivation", "Sample_Group"),
                   )

genes <- unique(c(all.features$"exp1_XGBoost",
                  all.features$"exp2_XGBoost",
                  all.features$"exp3_XGBoost"))


XG.1 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A2$Sample_Name],
                   annotations = annotations_A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames=T,
                   colAnno = c("Supposed biallelic inactivation", "Sample_Group"),
                   )

XG.2 <- preComputePheatmap(data=data_Promoters[genes,
                                       annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   showRownames=T,
                   colAnno = c("Supposed biallelic inactivation", "Sample_Group"),
                   )



data.frame(cutree(XG.1$tree_col,2))["M3797110001", ]

annotations[annotations$Sample_Name == "M3797110001",]

plot.sel.features <- patchwork::wrap_plots(list(as.ggplot(LR.1), as.ggplot(LR.2),
                                                as.ggplot(RF.1), as.ggplot(RF.2),
                                                as.ggplot(XG.1), as.ggplot(XG.2)),
                                           ncol = 2) + plot_layout(guides = "collect")

ggsave(filename = "~/ATM_Analysis/svn/analyse/results/FinalMethylation/Classifier/final_analyses/heatmaps_features_selected.png",
       plot = plot.sel.features,
       dpi=300,
       scale=1,
       units="in",
       width = 20,
       height = 15)




from_list <- function(dataframe) {
    members = unique(unlist(as.list(dataframe), use.names=F))
    members = members[!is.na(members)] 
    return(data.frame(sapply(dataframe, function(set) members %in% set),
                      row.names = members))
}


matrix = from_list(all.features)
matrix$gene_name = rownames(matrix)
head(matrix)


library(ComplexUpset)
upset(
    matrix,
    intersect=colnames(matrix)[1:length(colnames(matrix)) - 1],
    base_annotations=list(
        'Intersection size'=(
            intersection_size(
                bar_number_threshold=1,
                color='grey9',
                fill='grey80'
            )
            + geom_text(
                  mapping=aes(label=gene_name),
                  size = 2.8,
                  position=position_stack(),
                  na.rm=TRUE,
                  vjust=2.5
              )
        )
    ),
    width_ratio=0.15,
    height_ratio=1/4
)



classif <- read_xlsx(file.path(Origin, '/Classifier/pathways/A2_Pathways_evaluation_feature_selection.xlsx'))

classif <- fill_merged(classif, "FS method")
classif <- fill_merged(classif, "Classifier")
classif <- classif[c("FS method", "Classifier", "test_Matthews")]
colnames(classif)[colnames(classif) == "FS method"] <- "pathway"
classif <- classif %>%
    group_by(pathway, Classifier) %>%
    summarise(mean = mean(test_Matthews),
              std = sd(test_Matthews))

plot <- ggplot(
    classif[classif$Classifier %in% c("Random Forest",
                                      "XGBoost", "Support Vector Machine",
                                      "Logistic Regression"), ],
    aes(x = pathway, y = mean)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    geom_errorbar(aes(x = pathway,
                      ymin = mean - std,
                      ymax = mean + std)) +
    geom_hline(yintercept = 0.80, linetype=2, color = "red") + 
    ggtitle("Distribution of Matthews scores across Cross-Validation Folds") +
    xlab("Pathway") + ylab("Matthews score") + theme_bw() +
    theme(text=element_text(size=20),
          #axis.text.x = element_text(size=15),#angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) + 
    facet_wrap(~Classifier, as.table = FALSE, nrow=1) + coord_flip()

ggsave(filename = file.path(Origin, "/Classifier/pathways/A2_MCCS_pathways.png"), plot,
       width = 60,
       height = 50,
       units = "cm",
       dpi=300)


classif <- read_xlsx(file.path(Origin, '/Classifier/pathways/A2Inactive_Pathways_evaluation_feature_selection.xlsx'))

classif <- fill_merged(classif, "FS method")
classif <- fill_merged(classif, "Classifier")
classif <- classif[c("FS method", "Classifier", "test_Matthews")]
colnames(classif)[colnames(classif) == "FS method"] <- "pathway"
classif <- classif %>%
    group_by(pathway, Classifier) %>%
    summarise(mean = mean(test_Matthews),
              std = sd(test_Matthews))

plot <- ggplot(
    classif[classif$Classifier %in% c("Random Forest",
                                      "XGBoost", "Support Vector Machine",
                                      "Logistic Regression"), ],
    aes(x = pathway, y = mean)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    geom_errorbar(aes(x = pathway,
                      ymin = mean - std,
                      ymax = mean + std)) +
    geom_hline(yintercept = 0.80, linetype=2, color = "red") + 
    ggtitle("Distribution of Matthews scores across Cross-Validation Folds") +
    xlab("Pathway") + ylab("Matthews score") + theme_bw() +
    theme(text=element_text(size=20),
          #axis.text.x = element_text(size=15),#angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) + 
    facet_wrap(~Classifier, as.table = FALSE, nrow=1) + coord_flip()

ggsave(filename = file.path(Origin, "/Classifier/pathways/A2Inactive_MCCS_pathways.png"), plot,
       width = 60,
       height = 50,
       units = "cm",
       dpi=300)



# All promoters
classif <- read_xlsx(file.path(Origin, '/Classifier/A2_evaluation_several_ML.xlsx'))

classif <- fill_merged(classif, "Classifier")
classif <- classif[c("Classifier", "Fold", "test_Matthews")]


classif <- classif %>%
    group_by(Classifier) %>%
    summarise(mean = mean(test_Matthews),
              std = sd(test_Matthews))

head(classif)

plot <- ggplot(
    classif[classif$Classifier %in% c("Random Forest",
                                      "XGBoost", "Support Vector Machine",
                                      "Logistic Regression"), ],
    aes(x = Classifier, y = mean, fill = Classifier)) +
    geom_bar(stat = "identity", ) +
    geom_errorbar(aes(x = Classifier,
                      ymin = mean - std,
                      ymax = mean + std)) +
    geom_hline(yintercept = 0.80, linetype=2, color = "red") + 
    ggtitle("Distribution of Matthews scores across Cross-Validation Folds") +
    xlab("Classifier") + ylab("Matthews score") + theme_bw() +
    theme(text=element_text(size=20),
          #axis.text.x = element_text(size=15),#angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) #+ 
#    facet_wrap(~Classifier, as.table = FALSE, nrow=1) + coord_flip()

ggsave(filename = file.path(Origin, "/Classifier/A2_MCCS_genome-wide.png"), plot,
       width = 60,
       height = 50,
       units = "cm",
       dpi=300)







classif <- read_xlsx(file.path(Origin, '/Classifier/A2Inactive_evaluation_several_ML.xlsx'))

classif <- fill_merged(classif, "Classifier")
classif <- classif[c("Classifier", "Fold", "test_Matthews")]


classif <- classif %>%
    group_by(Classifier) %>%
    summarise(mean = mean(test_Matthews),
              std = sd(test_Matthews))

head(classif)

plot <- ggplot(
    classif[classif$Classifier %in% c("Random Forest",
                                      "XGBoost", "Support Vector Machine",
                                      "Logistic Regression"), ],
    aes(x = Classifier, y = mean, fill = Classifier)) +
    geom_bar(stat = "identity", ) +
    geom_errorbar(aes(x = Classifier,
                      ymin = mean - std,
                      ymax = mean + std)) +
    geom_hline(yintercept = 0.80, linetype=2, color = "red") + 
    ggtitle("Distribution of Matthews scores across Cross-Validation Folds") +
    xlab("Classifier") + ylab("Matthews score") + theme_bw() +
    theme(text=element_text(size=20),
          #axis.text.x = element_text(size=15),#angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) #+ 
#    facet_wrap(~Classifier, as.table = FALSE, nrow=1) + coord_flip()

ggsave(filename = file.path(Origin, "/Classifier/A2Inactive_MCCS_genome-wide.png"), plot,
       width = 60,
       height = 50,
       units = "cm",
       dpi=300)
