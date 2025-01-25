#!/usr/bin/env Rscript
#Torque Configuration
#PBS -l walltime=20:10:30
#PBS -l mem=100gb
#PBS -l nodes=2:ppn=4
#PBS -q batch
#PBS -N Mapping
#PBS -j oe

print("beginning of the script")

#---- Script du 25.10.2021 ----#
#---- Script minimalisé pour réaliser différentes analyses différentielles après avoir chargé les donnnées de méthylation ----#

#---- Script du 18.09.2023 ----#
#---- Script minimalisé afin de séparer le mapping des analyses différentielles  ----#

source("~/ATM_Analysis/svn/analyse/script/methylation/MethPipeline/config.R")

# Activate renv project
renv::activate(project = "~/ATM_Analysis/svn/analyse/script/")

# Lien vers les annotations et les Beta values
#Origin <- "/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation/"

#library("optparse")
#option_list = list(
#    make_option(c("-d", "--outDir"), type="character", default=NULL,
#                help="output directory", metavar="character")
#)
#opt_parser = OptionParser(option_list=option_list)
#opt = parse_args(opt_parser)

#Origin = opt$outDir

#Origin = "/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation/"
DirFiles = file.path(Origin, "matrices/")

# Liens vers les fichiers de sortie
#if(!dir.exists(file.path(Origin, "exploratoryAnalysis")))
#{
#    dir.create(file.path(Origin, "exploratoryAnalysis"), recursive = TRUE)
#}


# Lecture du fichier des M values (Les M values ayant été calculées dans le script précédent et les valeurs infinies supprimées)
load(file.path(DirFiles, "M_Funnorm_Normalisation.RData"))
load(file.path(DirFiles, "Beta_Funnorm_Normalisation.RData"))

# chargement des annotations
annotations <- read.csv(file.path(DirFiles, "Annotations.csv"), row.names = 1)

# Chargement des fonctions customisées
if (TCGA == FALSE){
    source(file = file.path(function.path, "functions.R"))
}

library("biomaRt")
library("ggplot2")

# Chargement des annotations des sondes de la puce EPIC
#library(ENmix)
#mf="~/ATM_Analysis/data/methylation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"
#annEPIC <- readmanifest(mf)
#annEPIC <- annEPIC$assay


GetFinalMatrix <- function(type, which_MSet, which_Feature, correspondance_table){
    require("dplyr")
    #Récupérer les M values ou les Beta value
    if (type == "Beta"){
        matrice  <- data.frame(all_beta)
    }else{
        matrice <- data.frame(all_M)
    }
    if (is.data.frame(correspondance_table) == TRUE){
        #Remplacer les Nan and Inf values by NA
        matrice[matrice == "-Inf" | matrice =="Inf" | is.na(matrice)] <- NA
        matrice$probe <- rownames(matrice)
        matrice <- merge(matrice, correspondance_table, by = "probe")
        countProbe <- matrice %>%
            group_by(SYMBOL) %>%
            summarise(NbrProbe = n())
#       print(head(countProbe))
        all_M_Final <- matrice %>%
            group_by(SYMBOL) %>%
            summarise_if(is.numeric, mean)
        all_M_Final <- data.frame(all_M_Final)
        all_M_Final <- merge(all_M_Final, countProbe, by = "SYMBOL")
        print(head(all_M_Final))
        rownames(all_M_Final) <- all_M_Final$SYMBOL
        all_M_Final <- subset(all_M_Final, select = -SYMBOL)
        return(all_M_Final)
    }else{
        matrice <-  subset(matrice, rownames(matrice) %in% correspondance_table)
        return(matrice)
    }
}


retry_with_delay <- function(expr, max_attempts = 50, delay_seconds = 30) {
  attempt <- 1
  while (attempt <= max_attempts) {
    try_result <- tryCatch(expr(), error = function(e) e)
    if (inherits(try_result, "error")) {
      cat("Attempt", attempt, "failed:", conditionMessage(try_result), "\n")
      cat("Retrying in", delay_seconds, "seconds...\n")
      Sys.sleep(delay_seconds)
      attempt <- attempt + 1
    } else {
      return(try_result)
    }
  }
  stop("Max attempts reached without success.")
}


connection <- function(){
# www, uswest, useast, asia
    ensembl <-  useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", mirror="www")
    return(ensembl)
}

ensembl <- retry_with_delay(connection)

mart <- useDataset(dataset="hsapiens_gene_ensembl", mart=ensembl)


# Number of probes present in the manifest
PROBES <- read.csv(file = "~/ATM_Analysis/data/annotation/updateFinal/Probes_annotations.bed", sep = "\t", header=F)
ProbesInManifest <- unique(PROBES[, "V4"])
NProbesInManifest <- length(ProbesInManifest)

# Number of probes present in the matrix of M values
NProbesMvalues <- length(rownames(all_M))


print("Genes")

# Chargement des sondes présentes dans les gènes
probes_Genes <- read.table(text = system(command="/bioinfo/local/build/Centos/bedtools/bedtools-2.25.0/bin/intersectBed -a ~/ATM_Analysis/data/annotation/updateFinal/Probes_annotations.bed -b ~/ATM_Analysis/data/annotation/updateFinal/Genes_annotations.bed -wa -wb", intern=TRUE))
probes_Genes <- probes_Genes[, c("V4", "V9")]
names(probes_Genes) <- c("probe", "ENSEMBL")
ProbesInGenes <- probes_Genes[, "probe"]
NGenesCorr <- dim(probes_Genes)[1]
NGenesCorrProbeUnique <- length(unique(probes_Genes$probe))
NGenesCorrENSEMBLUnique <-  length(unique(probes_Genes$ENSEMBL))


print(head(probes_Genes))

# Mapp EMSEMBL to SYMBOL
toSymbol <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=probes_Genes$ENSEMBL,mart= mart)
names(toSymbol) <-  c("ENSEMBL", "SYMBOL")
toSymbol <- subset(toSymbol, ! toSymbol$SYMBOL == "")
NGenesMapped <- length(unique(toSymbol$ENSEMBL))
PercentGenesMapped <- scales::percent(NGenesMapped / NGenesCorrENSEMBLUnique, accuracy = 0.01)

# Merge annotations and select only probe/Symbol correspondance
probes_Genes <- merge(x = probes_Genes,
      y = toSymbol,
      by = "ENSEMBL")
probes_Genes <- probes_Genes[, c('probe', 'SYMBOL')]

NProbesGenesFinal <- length(unique(probes_Genes$probe))
NGenesFinal <- length(unique(probes_Genes$SYMBOL))


print("Promoters")

# Chargement des sondes présentes dans les promoteurs
probes_Promoters <- read.table(text = system(command="/bioinfo/local/build/Centos/bedtools/bedtools-2.25.0/bin/intersectBed -a ~/ATM_Analysis/data/annotation/updateFinal/Probes_annotations.bed -b ~/ATM_Analysis/data/annotation/updateFinal/Promoters_annotations.bed -wa -wb", intern=TRUE))
probes_Promoters <- probes_Promoters[, c("V4", "V9")]
names(probes_Promoters) <- c("probe", "ENSEMBL")
ProbesInPromoters <- probes_Promoters[, "probe"]
NPromotersCorr <-  dim(probes_Promoters)[1]
NPromotersCorrProbeUnique <-  length(unique(probes_Promoters$probe))
NPromotersCorrENSEMBLUnique <-  length(unique(probes_Promoters$ENSEMBL))

# Mapp EMSEMBL to SYMBOL
toSymbolPromoters <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=probes_Promoters$ENSEMBL,mart= mart)
names(toSymbolPromoters) <-  c("ENSEMBL", "SYMBOL")
toSymbolPromoters <- subset(toSymbolPromoters, ! toSymbolPromoters$SYMBOL == "")

NPromotersMapped <- length(unique(toSymbolPromoters$ENSEMBL))
PercentPromotersMapped <- scales::percent(NPromotersMapped / NPromotersCorrENSEMBLUnique, accuracy = 0.01)

# Merge annotations and select only probe/Symbol correspondance
probes_Promoters <- merge(x = probes_Promoters,
                          y = toSymbolPromoters,
                          by = "ENSEMBL")
probes_Promoters <- probes_Promoters[, c('probe', 'SYMBOL')]

NProbesPromotersFinal <- length(unique(probes_Promoters$probe))
NPromotersFinal <- length(unique(probes_Promoters$SYMBOL))


print("Genes and promoters")

# Chargement des sondes présentes dans les gènes et promoteurs
probes_GenesAndPro <- read.table(text = system(command="/bioinfo/local/build/Centos/bedtools/bedtools-2.25.0/bin/intersectBed -a ~/ATM_Analysis/data/annotation/updateFinal/Probes_annotations.bed -b ~/ATM_Analysis/data/annotation/updateFinal/GeneAndPromoters_annotations.bed -wa -wb", intern=TRUE))
probes_GenesAndPro <- probes_GenesAndPro[, c("V4", "V9")]
names(probes_GenesAndPro) <- c("probe", "ENSEMBL")
ProbesInGenesAndPro <- probes_GenesAndPro[, "probe"]
NGenesAndProCorr <-  dim(probes_GenesAndPro)[1]
NGenesAndProCorrProbeUnique <-  length(unique(probes_GenesAndPro$probe))
NGenesAndProCorrENSEMBLUnique <-  length(unique(probes_GenesAndPro$ENSEMBL))

# Mapp EMSEMBL to SYMBOL
toSymbolGenesAndPro <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"),values = probes_GenesAndPro$ENSEMBL, mart = mart)
names(toSymbolGenesAndPro) <-  c("ENSEMBL", "SYMBOL")
toSymbolGenesAndPro <- subset(toSymbolGenesAndPro, ! toSymbolGenesAndPro$SYMBOL == "")

NGenesAndProMapped <- length(unique(toSymbolGenesAndPro$ENSEMBL))
PercentGenesAndProMapped <- scales::percent(NGenesAndProMapped / NGenesAndProCorrENSEMBLUnique, accuracy = 0.01)

# Merge annotations and select only probe/Symbol correspondance
probes_GenesAndPro <- merge(x = probes_GenesAndPro,
                            y = toSymbolGenesAndPro,
                            by = "ENSEMBL")
probes_GenesAndPro <- probes_GenesAndPro[, c('probe', 'SYMBOL')]

NProbesGenesAndProFinal <- length(unique(probes_GenesAndPro$probe))
NGenesAndProFinal <- length(unique(probes_GenesAndPro$SYMBOL))


print("CpG islands")

# Chargement des sondes présentes dans les CpG Islands
CpGIslands_annotations <- read.csv("~/ATM_Analysis/data/annotation/updateFinal/CpGIsland_annotations.bed", sep="\t", header=F)
CpGIslands_annotations$V4 <- with(CpGIslands_annotations, ave(as.character(V4), FUN = make.unique))
names(CpGIslands_annotations) <- c("chr", "start", "end", "CpG")
# Modification et sauvegarde
write.table(CpGIslands_annotations,
          "~/ATM_Analysis/data/annotation/updateFinal/CpGIsland_annotations_modified.bed", row.names=F, col.names=F, sep="\t", quote=F)

probes_CpGIsland <- read.table(text = system(command="/bioinfo/local/build/Centos/bedtools/bedtools-2.25.0/bin/intersectBed -a ~/ATM_Analysis/data/annotation/updateFinal/Probes_annotations.bed -b /data/users/nviart/ATM_Analysis/data/annotation/updateFinal/CpGIsland_annotations_modified.bed -wa -wb", intern=TRUE))
probes_CpGIsland <- probes_CpGIsland[, c("V4", "V9")]
names(probes_CpGIsland) <- c("probe", "SYMBOL")
ProbesInCpGIsland <- probes_CpGIsland[, "probe"]


print("Repeats")

# Chargement des sondes présentes dans les repeats
repeats_annotations <- read.csv("~/ATM_Analysis/data/annotation/update/Repeats_annotations.bed", sep="\t", header=F)
repeats_annotations$V5 <- with(repeats_annotations, ave(as.character(V5), FUN = make.unique))
names(repeats_annotations) <- c("chr", "start", "end", "strand", "CpG")
# Modification et sauvegarde
write.table(repeats_annotations,
          "~/ATM_Analysis/data/annotation/updateFinal/repeats_annotations_modified.bed", row.names=F, col.names=F, sep="\t", quote=F)

probes_Repeats <- read.table(text = system(command="/bioinfo/local/build/Centos/bedtools/bedtools-2.25.0/bin/intersectBed -a ~/ATM_Analysis/data/annotation/updateFinal/Probes_annotations.bed -b /data/users/nviart/ATM_Analysis/data/annotation/updateFinal/repeats_annotations_modified.bed -wa -wb", intern=TRUE))
probes_Repeats <- probes_Repeats[, c("V4", "V9")]
names(probes_Repeats) <- c("probe", "SYMBOL")
ProbesInRepeats <- probes_Repeats[, "probe"]




# Write all data
write.csv(probes_Genes, file.path(DirFiles, "probes_Genes.csv"))

write.csv(probes_Promoters, file.path(DirFiles, "probes_Promoters.csv"))

write.csv(probes_GenesAndPro, file.path(DirFiles, "probes_GenesAndPromoters.csv"))

write.csv(probes_Repeats, file.path(DirFiles, "probes_GenesAndPromoters.csv"))

write.csv(probes_CpGIsland, file.path(DirFiles, "probes_GenesAndPromoters.csv"))
 

plotsContext <- function(probeList,
                         annotations = NULL,
                         probes      = TRUE,
                         whichOne    = "first", # first, second, both
                         main        = "")
{
    require(ggplot2)
    require(patchwork)
    require(ggpubr)
    require(dplyr)

    if (is.null(annotations))
    {
        warning("You are using all samples of the annotations")
        data <- all_M
        anno <- annotations
    } else {
        data <- all_M[, annotations[["Sample_Name"]]]
        anno <- annotations
    }
    # subset probes
    if (probes)
    {
        data <- data[rownames(data) %in% probeList, ]
    } else {
        data <- probeList
    }

#    anno <- eval(parse(text=anno))
    
    # number of probes
    nprobes <- length(rownames(data))

    #print(head(data))

    data <- t(data)
    data <- data.frame(row.names = rownames(data), mean = rowMeans(data))
    if ("Variant_type_aggregated" %in% colnames(anno)){
        data <- merge(data, anno[c("Sample_Group", "Variant_type_aggregated")],
                      by="row.names")
    } else {
        data <- merge(data, anno[c("Sample_Group")],
                      by="row.names")
    }
    data[data$Sample_Group == "Non.ATM", ]$Sample_Group <- "Sporadic"
    data$Sample_Group <- factor(data$Sample_Group, levels = c("ATM", "Sporadic"))
    max <- max(data$mean) + 0.2

    require(rstatix)
    print(data %>% group_by(Sample_Group) %>%
                shapiro_test(mean))
    # custom theme
    custom.theme <- theme(legend.position="none",
                          panel.grid = element_line(colour = "grey92"),
                          panel.grid.minor = element_line(linewidth = rel(0.5)),
                          panel.grid.major = element_line(linewidth = rel(0.5)))

    col <- c("ATM" = "#F8756E", "Sporadic" = "#01BEC5")
    if (whichOne %in% c("first", "both"))
    {
    # First comparison
#        wilcox.test (…paired = FALSE) /// kruskal.test /// oneway <- test       
        comparisons <- list( c("Sporadic", "ATM") )
        t.test <- compare_means(mean ~ Sample_Group, comparisons = comparisons,
                                p.adj= "BH", method='wilcox.test', data = data)
        t.text <- t.test %>%
            mutate(y.position = c(-1.5))
        t.test$label <- paste0(t.test$p.adj, " (", t.test$p.signif, ")")
        print(t.test)

        data <- data %>% group_by(Sample_Group) %>% mutate(count=n()) %>% data.frame()
        #data$Sample_Group <- paste0(data$Sample_Group, "\n(", data$count, ")")
        #print(data)
        
        p1 <- ggviolin(data,
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
    }

#    if (whichOne %in% c("second", "both"))
#    {
#    # Second comparison
#        comparisons <- list( c("LoF", "MV"), c("MV", "Non.carrier"), c("LoF", "Non.carrier") )
#        data$Variant_type_aggregated[data$Variant_type_aggregated == ""] <- "Non carrier"
#        
#        t.test <- compare_means(mean ~ Variant_type_aggregated,
#                                comparisons = comparisons,
#                                p.adj= "BH", method='t.test', data = data)
#        t.test <- t.test %>%
#            mutate(y.position = c(-1.5))
#        t.test$label <- paste0(t.test$p.adj, " (", t.test$p.signif, ")")
#        
#        p2 <- ggviolin(data,
#                       x = "Variant_type_aggregated",
#                       y = "mean",
#                       palette = c("#01BEC5", "#FF9D98", "#B8544F"),
#                       fill = "Variant_type_aggregated",
#                       add = "boxplot") +
#            theme(legend.position="none") +
#            stat_pvalue_manual(t.test, label = "label",
#                               y.position = c(max + 0.2, max, max + 0.2)) +
#            xlab("Type of variant") +
#            custom.theme
#    }

    if (whichOne == "first")
    {
        return(p1)
    } else {
        if (whichOne == "second"){
            return(p2)
        } else {
            allPlots <- ggarrange(p2, p2 + theme_bw(),
                                  labels = c("A", "B"),
                                  ncol = 2)
            allPlots <- annotate_figure(allPlots,
                                        top = text_grob(main, size=18))
            return(allPlots)
        }
    }    
}


comment <- function(){

ProbesInPromoters <- probes_Promoters[, 1]
ProbesInGenes <- probes_Genes[, 1]
ProbesInCpGIsland <- probes_CpGIsland[, 1]
ProbesInRepeats <- probes_Repeats[, 1]
    
p1.all <- plotsContext(rownames(all_M),
                       whichOne = "first",
                       anno = annotations_A0A2,
                       main = "All probes")

p1.pro <- plotsContext(ProbesInPromoters,
                       whichOne = "first",
                       anno = annotations_A0A2,
                       main = "Probes mapped to promoters")

p1.gene <- plotsContext(ProbesInGenes,
                       whichOne = "first",
                       anno = annotations_A0A2,
                       main = "Probes mapped to genes")

p1.cpg <- plotsContext(ProbesInCpGIsland,
                       whichOne = "first",
                       anno = annotations_A0A2,
                       main = "Probes mapped to CpG islands")

p1.repeats <- plotsContext(ProbesInRepeats,
                           whichOne = "first",
                           anno = annotations_A0A2,
                           main = "Probes mapped to repeats")


allPlots <- ggarrange(p1.all, p1.pro, p1.gene, p1.cpg, p1.repeats,
                      labels = c("A", "B", "C", "D", "E"),
                      ncol = 3, nrow = 2) + bgcolor("White")   
allPlots <- annotate_figure(allPlots,
                            top = text_grob("Mean level of methylation according to the tumour type for:", size=18))

ggsave(filename = file.path(Origin, "exploratoryAnalysis/A0A2_genomic_context.png"), allPlots, width=18, height=14, bg="white")

       

p1.all <- plotsContext(rownames(all_M),
                       whichOne = "first",
                       anno = annotations_A2,
                       main = "All probes")

p1.pro <- plotsContext(ProbesInPromoters,
                       whichOne = "first",
                       anno = annotations_A2,
                       main = "Probes mapped to promoters")

p1.gene <- plotsContext(ProbesInGenes,
                       whichOne = "first",
                       anno = annotations_A2,
                       main = "Probes mapped to genes")

p1.cpg <- plotsContext(ProbesInCpGIsland,
                       whichOne = "first",
                       anno = annotations_A2,
                       main = "Probes mapped to CpG islands")

p1.repeats <- plotsContext(ProbesInRepeats,
                           whichOne = "first",
                           anno = annotations_A2,
                           main = "Probes mapped to repeats")


allPlots <- ggarrange(p1.all, p1.pro, p1.gene, p1.cpg, p1.repeats,
                      labels = c("A", "B", "C", "D", "E"),
                      ncol = 3, nrow = 2) + bgcolor("White")   
allPlots <- annotate_figure(allPlots,
                            top = text_grob("Mean level of methylation according to the tumour type for:", size=18))

ggsave(filename = file.path(Origin, "exploratoryAnalysis/A2_genomic_context.png"), allPlots, width=18, height=14, bg="white")






p1.all <- plotsContext(rownames(all_M),
                       whichOne = "first",
                       anno = annotations_AVUSA2,
                       main = "All probes")

p1.pro <- plotsContext(ProbesInPromoters,
                       whichOne = "first",
                       anno = annotations_AVUSA2,
                       main = "Probes mapped to promoters")

p1.gene <- plotsContext(ProbesInGenes,
                       whichOne = "first",
                       anno = annotations_AVUSA2,
                       main = "Probes mapped to genes")

p1.cpg <- plotsContext(ProbesInCpGIsland,
                       whichOne = "first",
                       anno = annotations_AVUSA2,
                       main = "Probes mapped to CpG islands")

p1.repeats <- plotsContext(ProbesInRepeats,
                           whichOne = "first",
                           anno = annotations_AVUSA2,
                           main = "Probes mapped to repeats")


allPlots <- ggarrange(p1.all, p1.pro, p1.gene, p1.cpg, p1.repeats,
                      labels = c("A", "B", "C", "D", "E"),
                      ncol = 3, nrow = 2) + bgcolor("White")   
allPlots <- annotate_figure(allPlots,
                            top = text_grob("Mean level of methylation according to the tumour type for:", size=18))

ggsave(filename = file.path(Origin, "exploratoryAnalysis/AVUSA2_genomic_context.png"), allPlots, width=18, height=14, bg="white")




}






# Get final matrices with the number of probes


# Get the final matrices with the number of probes
data_Genes <- GetFinalMatrix("M", all_M, "genes", probes_Genes)
data_Promoters <- GetFinalMatrix("M", all_M, "promoters", probes_Promoters)
data_GenesAndPro <- GetFinalMatrix("M", all_M, "GenesAndPro", probes_GenesAndPro)

data_CpGIsland <- GetFinalMatrix("M", all_M, "CpGIsland", probes_CpGIsland)

data_Repeats <- GetFinalMatrix("M", all_M, "Repeats", probes_Repeats)

# Write the number of probes summarised
write.csv(data_Genes["NbrProbe"], file.path(DirFiles, "M_Genes_NbrProbes.csv"))
write.csv(data_Promoters["NbrProbe"], file.path(DirFiles, "M_Promoterss_NbrProbes.csv"))
write.csv(data_GenesAndPro["NbrProbe"], file.path(DirFiles, "M_GenesAndPromoters_NbrProbes.csv"))





ProbesInM <- unique(rownames(all_M))

NbrProbeMGenes <- length(unique(probes_Genes[probes_Genes$probe %in% ProbesInM, ]$probe))
NGenesMFinal <- length(unique(probes_Genes[probes_Genes$probe %in% ProbesInM, ]$SYMBOL))

NbrProbeMPromoters <- length(unique(probes_Promoters[probes_Promoters$probe %in% ProbesInM, ]$probe))
NPromotersMFinal <- length(unique(probes_Promoters[probes_Promoters$probe %in% ProbesInM, ]$SYMBOL))

NbrProbeMGenesAndPro <- length(unique(probes_GenesAndPro[probes_GenesAndPro$probe %in% ProbesInM, ]$probe))
NGenesAndProMFinal <- length(unique(probes_GenesAndPro[probes_GenesAndPro$probe %in% ProbesInM, ]$SYMBOL))




# Suppress the number of probes
data_Genes <- subset(data_Genes, select = -NbrProbe)
data_Promoters <- subset(data_Promoters, select = -NbrProbe)
data_GenesAndPro <- subset(data_GenesAndPro, select = -NbrProbe)

data_CpGIsland <- subset(data_CpGIsland, select = -NbrProbe)

data_Repeats <- subset(data_Repeats, select = -NbrProbe)

# Write data
write.csv(data_Genes, file.path(DirFiles, "M_Genes_annotated.csv"))
write.csv(data_Promoters, file.path(DirFiles, "M_Promoters_annotated.csv"))
write.csv(data_GenesAndPro, file.path(DirFiles, "M_GenesAndPromoters_annotated.csv"))

write.csv(data_CpGIsland, file.path(DirFiles, "M_CpGIsland_annotated.csv"))

write.csv(data_Repeats, file.path(DirFiles, "M_Repeats_annotated.csv"))


# Get the final matrices with the number of probes
data_Genes_Beta <- GetFinalMatrix("Beta", all_Beta, "genes", probes_Genes)
data_Promoters_Beta <- GetFinalMatrix("Beta", all_Beta, "promoters", probes_Promoters)
data_GenesAndPro_Beta <- GetFinalMatrix("Beta", all_Beta, "GenesAndPro", probes_GenesAndPro)

# Write the number of probes summarised
write.csv(data_Genes_Beta["NbrProbe"], file.path(DirFiles, "Beta_Genes_NbrProbes.csv"))
write.csv(data_Promoters_Beta["NbrProbe"], file.path(DirFiles, "Beta_Promoterss_NbrProbes.csv"))
write.csv(data_GenesAndPro_Beta["NbrProbe"], file.path(DirFiles, "Beta_GenesAndPromoters_NbrProbes.csv"))

# Suppress the number of probes
data_Genes_Beta <- subset(data_Genes_Beta, select = -NbrProbe)
data_Promoters_Beta <- subset(data_Promoters_Beta, select = -NbrProbe)
data_GenesAndPro_Beta <- subset(data_GenesAndPro_Beta, select = -NbrProbe)

# Write data
write.csv(data_Genes_Beta, file.path(DirFiles, "Beta_Genes_annotated.csv"))
write.csv(data_Promoters_Beta, file.path(DirFiles, "Beta_Promoters_annotated.csv"))
write.csv(data_GenesAndPro_Beta, file.path(DirFiles, "Beta_GenesAndPromoters_annotated.csv"))





save(probes_Genes, probes_Promoters, probes_GenesAndPro, probes_Repeats, probes_CpGIsland, data_Genes_Beta, data_Promoters_Beta, data_GenesAndPro_Beta, data_Genes, data_Promoters, data_GenesAndPro, data_CpGIsland, data_Repeats,
     file = file.path(Origin, "matrices/Mapping.RData"))





# Output the table describing all stats 
colNames <- c("Number of probes", "Number of targeted element", "percentage mapped (in %)", "Number of probes", "Number of targeted element", "Number of probes", "percentage of final probes (in %)", "Number of targeted element")
#scales::percent(NGenesCorrProbeUnique/NProbesMvalues, accuracy = 0.01)
GenesRow <- c(NGenesCorrProbeUnique, NGenesCorrENSEMBLUnique, PercentGenesMapped, NProbesGenesFinal, NGenesFinal, NbrProbeMGenes, scales::percent(NbrProbeMGenes/NProbesMvalues, accuracy = 0.01), NGenesMFinal)
#scales::percent(NPromotersCorrProbeUnique/NProbesMvalues, accuracy = 0.01)
PromotersRow <- c(NPromotersCorrProbeUnique, NPromotersCorrENSEMBLUnique, PercentPromotersMapped, NProbesPromotersFinal, NPromotersFinal, NbrProbeMPromoters, scales::percent(NbrProbeMPromoters/NProbesMvalues, accuracy = 0.01), NPromotersMFinal)
#scales::percent(NGenesAndProCorrProbeUnique/NProbesMvalues, accuracy = 0.01)
GenesAndProRow <- c(NGenesAndProCorrProbeUnique, NGenesAndProCorrENSEMBLUnique, PercentGenesAndProMapped, NProbesGenesAndProFinal, NGenesAndProFinal, NbrProbeMGenesAndPro, scales::percent(NbrProbeMGenesAndPro/NProbesMvalues, accuracy = 0.01), NGenesAndProMFinal)

dataFrame <- rbind(GenesRow, PromotersRow, GenesAndProRow)
colnames(dataFrame) <- colNames
rownames(dataFrame) <- c("Genes", "Promoters", "Genes and Promoters")




# Export the table
knitr::kable(dataFrame, "html",
             booktabs = TRUE,
             caption = "<center><strong>Summary of genome contexts covered by the Illumina Infinium EPIC BeadChip</strong></center>",
             escape = FALSE,
             align = c('l', rep('c', 11))) %>%
    kableExtra::kable_styling(latex_options = "striped", bootstrap_options = c("striped", "hover", "condensed", "responsive", fixed_thead = T), full_width = F) %>%
                        kableExtra::add_header_above(c(" " = 1,
                                                       setNames(2, "ENSEMBL"),
                                                       setNames(1, "Correspondance \t ENSEMBL/SYMBOL"),
                                                       setNames(2, "SYMBOL"),
                                                     setNames(3, "Final matrix of M values")))  %>%             
        kableExtra::save_kable(file.path(Origin, "Genomic_context_stat.html"))














comment <- function(){

# Chargement des sondes présentes dans les repeats 
probes_Repeats <- read.table(text = system(command="/bioinfo/local/build/Centos/bedtools/bedtools-2.25.0/bin/intersectBed -a ~/ATM_Analysis/data/annotation/updateFinal/Probes_annotations.bed -b /data/users/nviart/ATM_Analysis/data/annotation/update/Repeats_annotations.bed -wa -wb", intern=TRUE))
probes_Repeats <- probes_Repeats[, c("V4", "V9", "V10")]
names(probes_Repeats) <- c("probe", "CpG", "context")
ProbesInRepeats <- probes_Repeats[, "probe"]


res = list()
var.list <- unique(probes_Repeats$context)
for (i in 1:length(var.list)){
    print(var.list[i])
    d <- probes_Repeats[probes_Repeats$context == var.list[i], ]
    p <- plotsContext(d$probe,
                      whichOne = "first",
                      anno = annotations_A1,
                      main = paste0(var.list[i])
                      )
    res[[i]] <- p
}

allRepeats <- ggarrange(
    p1.repeats,
    res[[1]],
    res[[2]],
    res[[3]],
    res[[4]],
    res[[5]],
    res[[6]],
    res[[7]],
    res[[8]],
    res[[9]],
    res[[10]],
    res[[11]],
    res[[12]],
    res[[13]],
    res[[14]],
    res[[15]],
    res[[16]],
    res[[17]],
    res[[18]],
    res[[19]],
#    res[[20]],
    labels = LETTERS[1:21],
    ncol = 5, nrow = 4)


allRepeats

plotsContext(ProbesInPromoters,
             main = "Level of methylation of probes mapped to promoters")

plotsContext(ProbesInGenes)

plotsContext(ProbesInCpGIsland)

plotsContext(ProbesInRepeats)

plotsContext(PromotersCpG)

plotsContext(PromotersLowCpG)


#}

OutFiles = file.path(Origin, "differentialAnalysis/GenomeWide/")

names(probes_CpGIsland) <- c("probe", "SYMBOL")
data_CpGIslands <- GetFinalMatrix("M", all_M, "CpGs", probes_CpGIsland)
data_CpGIslands <- subset(data_CpGIslands, select = -NbrProbe)

#plotsContext(data_CpGIslands,
#             probes = FALSE,
#             whichOne = "first",
#             main = "")


# Differentiall analysis of CpG (A1 and A2) and retrieve thei position in the genome
result = list()
for (A in c("A1", "A2"))
{
    DifferentialAnalysis(data            = data_CpGIslands,
                         annotations     = eval(parse(text=paste0("annotations_", A))),
                         colInterest     = "Sample_Group",
                         analysisLevel   = "GenomeWide", #Gene or GenomeWide
                                        #listProbes      = ProbesInATMGene ,
                         contrastToUse   = "ATM - Non.ATM",
                         comparison_text = "Comparison of ATM and Non.ATM cancers",
                         FC.threshold    = 1,
                         p.val.threshold = 0.05,
                         path            = file.path(OutFiles,
                                                     paste0(A, "_CpGIslands_1FC_"))   ,
                         GO              = FALSE,
                         KEGG            = FALSE,
                         keepAll         = TRUE
                         )
    # read the CpG coordinates
    cpg <- read.csv("~/ATM_Analysis/data/annotation/updateFinal/CpGIsland_annotations_modified.bed", sep="\t", header=FALSE)
    names(cpg) <- c("chr_CpG", "start_CpG", "end_CpG", "CpG")
    # read the results of the differential analysis
    df <- read.csv(file = file.path(OutFiles,
                                    paste0(A, "_CpGIslands_1FC_DMPs.csv"))  ,
                   row.names = 1)
    # get only differentially methylated CpGs
    cpgs <- rownames(df[df$DifferentiallyMethylated == TRUE,])
    #names(cpg) <- c("chr", "start", "end", "CpG")
    annoCpG <- cpg[cpg$CpG %in% cpgs, ]
    # Add the coordinates of CpG
    annoCpG <- merge(df, annoCpG, by.x="row.names", by.y="CpG")
    names(annoCpG)[names(annoCpG) == 'Row.names'] <- 'CpG'
    annoCpG <- annoCpG[order(annoCpG$CpG), ]
    annoCpG$Analysis <- A
    result <- c(result, list(annoCpG))
}
result <- dplyr::bind_rows(result)

# write the result to file to intesect with genes
write.table(result[c("chr_CpG", "start_CpG", "end_CpG", "CpG", "Analysis", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val","B", "DifferentiallyMethylated")],
            file = file.path(OutFiles, "CpG_analysis.bed"),
            quote=F, row.names=F, col.names=F, sep="\t")

# Interect CpGIslands with gene positions
intersect <- read.table(text = system(command="/bioinfo/local/build/Centos/bedtools/bedtools-2.25.0/bin/intersectBed -a ~/ATM_Analysis/svn/analyse/results/FinalMethylation/differentialAnalysis/GenomeWide/CpG_analysis.bed -b ~/ATM_Analysis/data/annotation/updateFinal/Genes_annotations.bed -wa -wb", intern=TRUE))

names(intersect) <- c("chr_CpG", "start_CpG", "end_CpG", "CpG", "analysis", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val","B", "DifferentiallyMethylated", "chr", "start", "end", "ENSEMBL", "strand")

# Convert ENSEMBL to gene symbols
ensembl <-  useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", mirror="asia")
mart <- useDataset(dataset="hsapiens_gene_ensembl", mart=ensembl)

toSymbol <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=intersect$ENSEMBL, mart= mart)
names(toSymbol) <-  c("ENSEMBL", "SYMBOL")
toSymbol[toSymbol$SYMBOL == "", ]$SYMBOL <- toSymbol[toSymbol$SYMBOL == "", ]$ENSEMBL

# merge the results
result <- merge(intersect, toSymbol, by = "ENSEMBL")
result <- result[order(result$analysis, result$CpG), ]

# get annotation with the python script (3.2.annotate....)
write.csv(result, file.path(Origin, "differentialAnalysis/GenomeWide/CpG_analysis_genesFound.bed"))


#heatmap_func(data = data_CpGIslands[cpgs, annotations_A1$Sample_Name],
#             annotations       = annotations_A1,
#             subsetAnnotation  = c("Study", "ER_status", "Subtype", "Grade", "Variant_type", "ATM_LOH"),
#             color             = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
#             fontsize_row      = 8,
#             fontsize_col      = 10,
#             main              = 'test pheatmap',
#             )
}





