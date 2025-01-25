#!/usr/bin/env Rscript
#Torque Configuration
#PBS -l walltime=20:10:30
#PBS -l mem=60gb
#PBS -l nodes=2:ppn=4
#PBS -q batch
#PBS -N DifferentialAnalysis
#PBS -j oe

print("beginning of the script")

source("~/ATM_Analysis/svn/analyse/script/methylation/MethPipeline/config.R")

# Activate renv project
renv::activate(project = "~/ATM_Analysis/svn/analyse/script/")


#library("optparse")
#option_list = list(
#    make_option(c("-d", "--outDir"), type="character", default=NULL,
#                help="output directory", metavar="character")
#)
#opt_parser = OptionParser(option_list=option_list)
#opt = parse_args(opt_parser)

#Origin = file.path(opt$outDir)
#Origin = "/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation/"

# Lien vers les annotations et les Beta values
DirFiles = file.path(Origin, "matrices/")
# Liens vers les fihiers de sortie
OutFilesDA = file.path(Origin, "differentialAnalysis")

# Lecture du fichier des Beta values
#all_beta <- read.csv(file.path(DirFiles, "Beta_Funnorm_Normalisation.csv"), row.names = 1)
load(file.path(DirFiles, "Beta_Funnorm_Normalisation.RData"))

# Lecture du fichier des M values (Les M values ayant été calculées dans le script précédent et les valeurs infinies supprimées)
#all_M <- read.csv(file.path(DirFiles, "M_Funnorm_Normalisation.csv"), row.names = 1)
load(file.path(DirFiles, "M_Funnorm_Normalisation.RData"))
load(file.path(DirFiles, "Mapping.RData"))

# chargement des annotations
annotations <- read.csv(file.path(DirFiles, "Annotations.csv"), row.names = 1)

ATM.LOH.Yes <- c("T0220", "T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
ATM.LOH.No <- c("T0120", "T0015", "T0009", "T0249", "T0247")    
annotations[["ATM_LOH"]]  <- ifelse(annotations$Sample_Name %in% ATM.LOH.Yes, "LOH",
                             ifelse(annotations$Sample_Name %in% ATM.LOH.No, "No LOH",
                                    "Unknown LOH")
                             )


# Chargement des fonctions customisées
source(file = file.path(function.path, "functions.R"))


OutFilesAnno = file.path(Origin, "matrices/Annotations/")
if(!dir.exists(OutFilesAnno))
{
    dir.create(OutFilesAnno, recursive = TRUE)
}
write.csv(annotations_A0A2, file.path(OutFilesAnno, "annotations_A0A2.csv"))
write.csv(annotations_A2, file.path(OutFilesAnno, "annotations_A2.csv"))
write.csv(annotations_AVUSA2, file.path(OutFilesAnno, "annotations_AVUSA2.csv"))
write.csv(annotations_A2Inactive, file.path(OutFilesAnno, "annotations_A2Inactive.csv"))

# Chargement des annotations des sondes de la puce EPIC
#library(ENmix)
#mf="~/ATM_Analysis/data/methylation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"
#annEPIC <- readmanifest(mf)
#annEPIC <- annEPIC$assay

# Chargement des packages
#library(EnsDb.Hsapiens.v86)
library("biomaRt")
library("ggplot2")

library("scales")
library('bumphunter')
library("RColorBrewer")
library('limma')
library('DT')
library('plyr')
library('dplyr')
library('Biobase')
library('GEOquery')
library('data.table')
library("pheatmap")
library("gridExtra")
library("FactoMineR")
library("forcats")
library("backports")
library("factoextra")
library("ggpubr")
library("ggrepel")
library("clusterProfiler")
library("org.Hs.eg.db")
library("GenomeInfoDb")
library("RPMM")
library("lumi")
library("ENmix")
library("enrichplot")
library("cowplot")


# Stats on probes
p <- rownames(all_M)
p <- probes_Promoters[probes_Genes$probe %in% p, ]
print(cat("Number of genes: ", length(unique(p$SYMBOL))))
print(cat("Number of probes: ", length(unique(p$probe))))

p <- rownames(all_M)
p <- probes_Promoters[probes_Promoters$probe %in% p, ]
print(cat("Number of promoter: ", length(unique(p$SYMBOL))))
print(cat("Number of probes: ", length(unique(p$probe))))


test <- aggregate(SYMBOL ~., p, toString)

annoProbes <- read.csv("~/ATM_Analysis/data/annotation/updateFinal/Probes_annotations.bed", sep="\t", header=F)
colnames(annoProbes) <- c("chr", "start", "end", "probe", "strand")

testF <- merge(test, annoProbes, all.x=T, all.y=F)

write.csv(testF, file.path(Origin, "matrices/probes_Promoters_mapping.csv"))


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

ProbesInATMGene <- getProbesInGene("ATM", probes_Genes, "SYMBOL")
ProbesInATMPromoter <- getProbesInGene("ATM", probes_Promoters, "SYMBOL")

ProbesInBRCA1Gene <- getProbesInGene("BRCA1", probes_Genes, "SYMBOL")
ProbesInBRCA1Promoter <- getProbesInGene("BRCA1", probes_Promoters, "SYMBOL")
ProbesInBRCA2Gene <- getProbesInGene("BRCA2", probes_Genes, "SYMBOL")
ProbesInBRCA2Promoter <- getProbesInGene("BRCA2", probes_Promoters, "SYMBOL")





#########################################################################################
############################ Fonction Analyses différentielles ##########################
#########################################################################################

print("Beginning Differential Analysis on ATM")

######################
####  ATM analyses
######################

# Retrieve groups of annotations
#groups <- getGroups()
#annotations_A1 <- groups$annotations_A1
#annotations_A2 <- groups$annotations_A2
#annotations_A3 <- groups$annotations_A3
#annotations_A4 <- groups$annotations_A4
#annotations_Old <- groups$annotations_Old
#annotations_AOldER <- groups$annotations_AOldER
#annotations_AVUS <- groups$annotations_AVUS

OutFiles = file.path(OutFilesDA, "ByGene")

if(!dir.exists(OutFiles))
{
    dir.create(OutFiles, recursive = TRUE)
}


if (TCGA){
    rownames(annotations) <- stringr::str_replace_all(rownames(annotations), "-", ".")
}

# Analyses différentielles sur le gène et promoteur d'ATM
# Analyse A0: tous les échantillons
DifferentialAnalysis(data             = all_M,
                     annotations      = annotations,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMGene,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A0_ATM_Gene_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMPromoter ,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A0_ATM_Promoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = c(ProbesInATMGene, ProbesInATMPromoter)[!duplicated(c(ProbesInATMGene, ProbesInATMPromoter))],
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A0_ATM_GeneANDPromoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

# Analyse A1:  ATM pathogènes
DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A1,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMGene,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A1_ATM_Gene_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A1,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMPromoter ,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A1_ATM_Promoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A1,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = c(ProbesInATMGene, ProbesInATMPromoter)[!duplicated(c(ProbesInATMGene, ProbesInATMPromoter))],
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A1_ATM_GeneANDPromoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )



# Analyse AVUS:  ATM VUS
DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_AVUS,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMGene,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "AVUS_ATM_Gene_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_AVUS,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMPromoter ,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "AVUS_ATM_Promoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_AVUS,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = c(ProbesInATMGene, ProbesInATMPromoter)[!duplicated(c(ProbesInATMGene, ProbesInATMPromoter))],
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "AVUS_ATM_GeneANDPromoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )




## all ATM ER+
DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A0A2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMGene,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A0A2_ATM_Gene_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A0A2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMPromoter ,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A0A2_ATM_Promoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A0A2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = c(ProbesInATMGene, ProbesInATMPromoter)[!duplicated(c(ProbesInATMGene, ProbesInATMPromoter))],
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A0A2_ATM_GeneANDPromoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )





## all ATM PV ER+

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMGene,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A2_ATM_Gene_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMPromoter ,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A2_ATM_Promoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = c(ProbesInATMGene, ProbesInATMPromoter)[!duplicated(c(ProbesInATMGene, ProbesInATMPromoter))],
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A2_ATM_GeneANDPromoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

## all ATM VUS ER+


DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_AVUSA2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMGene,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "AVUSA2_ATM_Gene_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_AVUSA2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = ProbesInATMPromoter ,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "AVUSA2_ATM_Promoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_AVUSA2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "Gene", #Gene or GenomeWide
                     listProbes       = c(ProbesInATMGene, ProbesInATMPromoter)[!duplicated(c(ProbesInATMGene, ProbesInATMPromoter))],
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 0.5,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "AVUSA2_ATM_GeneANDPromoter_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )







print("Beginning of Genome Wide Differential Analysis")

######################
# Genome wide analyses
######################

OutFiles = file.path(OutFilesDA, "GenomeWide")
if(!dir.exists(OutFiles))
{
    dir.create(OutFiles, recursive = TRUE)
}


N.promoters <- read.csv(file = file.path(Origin, "/matrices/M_Promoterss_NbrProbes.csv"), row.names = 1)
#scale_values <- function(x){(x - min(x) + 0.1) / (max(x) - min(x))}
#weights.Pro <- scale_values(N.promoters[, 1])
weights.Pro <- N.promoters[, 1]

N.genes <- read.csv(file = file.path(Origin, "/matrices/M_Genes_NbrProbes.csv"), row.names = 1)
weights.Genes <- N.genes[, 1]


# Analyse A0: tous les échantillons
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Genes,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Pro,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )
# All probes
DifferentialAnalysis(data            = all_M,
                     annotations     = annotations,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0_Probes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


# Analyse A1: échantillons ATM pathogénique ou likely pathogéniques
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations_A1,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Genes,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A1_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations_A1,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Pro,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A1_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations_A1,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A1_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations_A1,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A1_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# All probes
DifferentialAnalysis(data            = all_M,
                     annotations     = annotations_A1,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A1_Probes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


# Analyse A2: échantillons ATM pathogénique ou likely pathogéniques
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations_A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Genes,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations_A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Pro,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations_A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations_A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# All probes
DifferentialAnalysis(data             = all_M,
                     annotations      = annotations_A2,
                     colInterest      = "Sample_Group",
                     analysisLevel    = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse    = "ATM - Non.ATM",
                     comparison_text  = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold     = 1,
                     p.val.threshold  = 0.05,
                     path             = file.path(OutFiles, "A2_Probes_"),
                     GO               = FALSE,
                     KEGG             = FALSE,
                     keepAll          = TRUE
                     )

comment <- function(){
# Analyse A3
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations_A3,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A3_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations_A3,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A3_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations_A3,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A3_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations_A3,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A3_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# All probes
DifferentialAnalysis(data            = all_M,
                     annotations     = annotations_A3,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A3_Probes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


# Analyse A4
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations_A4,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     ##listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A4_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations_A4,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     ##listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A4_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations_A4,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A4_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations_A4,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A4_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# All probes
DifferentialAnalysis(data            = all_M,
                     annotations     = annotations_A4,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     ##listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A4_Probes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )





# Analyse Old controls
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations_AOld,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOld_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations_AOld,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOld_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations_AOld,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOld_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations_AOld,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOld_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )
# All probes
DifferentialAnalysis(data            = all_M,
                     annotations     = annotations_AOld,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOld_Probes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )





# Analyse P or LP ATM variant & Old and ER+ controls
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations_AOldER,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOldER_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations_AOldER,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOldER_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations_AOldER,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOldER_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations_AOldER,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOldER_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )
# All probes
DifferentialAnalysis(data            = all_M,
                     annotations     = annotations_AOldER,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AOldER_Probes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

}








# Analyse AVUS: Variant of Uncertain Significance
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations_AVUS,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Genes,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUS_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations_AVUS,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Pro,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUS_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations_AVUS,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUS_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations_AVUS,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUS_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# All probes
DifferentialAnalysis(data            = all_M,
                     annotations     = annotations_AVUS,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUS_Probes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )












# Analyse AVUSA2: Variant of Uncertain Significance and ER+
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations_AVUSA2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Genes,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations_AVUSA2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Pro,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations_AVUSA2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations_AVUSA2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# All probes
DifferentialAnalysis(data            = all_M,
                     annotations     = annotations_AVUSA2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2_Probes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )






# Analyse A0A2: all ER+ tumours
# In Genes
DifferentialAnalysis(data            = data_Genes,
                     annotations     = annotations_A0A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Genes,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = annotations_A0A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Pro,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In CpGIslands
DifferentialAnalysis(data            = data_CpGIsland,
                     annotations     = annotations_A0A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2_CpGIsland_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Genes and Promoters
DifferentialAnalysis(data            = data_GenesAndPro,
                     annotations     = annotations_A0A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2_GenesAndPro_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# All probes
DifferentialAnalysis(data            = all_M,
                     annotations     = annotations_A0A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2_Probes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )




# Without samples with LOH

# In Promoters
DifferentialAnalysis(data            = data_Promoters,
                     annotations     = dplyr::filter(annotations_A0A2,
                                                     ATM_LOH %in% c("No LOH", "Unknown LOH")),
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Pro,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2_noLOH_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


DifferentialAnalysis(data            = data_Promoters,
                     annotations     = dplyr::filter(annotations_A2,
                                                     ATM_LOH %in% c("No LOH", "Unknown LOH")),
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Pro,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2_noLOH_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )




# With only samples with certain inactivation


# Analyse A2: échantillons ATM pathogénique ou likely pathogéniques
# In Genes
sel <- annotations_A2[(annotations_A2$Sample_Group == "ATM" &
                       annotations_A2$Inactive == "inactive ATM") |
                      annotations_A2$Sample_Group == "Non.ATM", ]$Sample_Name
        
DifferentialAnalysis(data            = data_Genes[sel],
                     annotations     = annotations_A2[sel, ],
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Genes,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2Inactive_Genes_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

# In Promoters
DifferentialAnalysis(data            = data_Promoters[sel],
                     annotations     = annotations_A2[sel, ],
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     weights         = weights.Pro,
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2Inactive_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )









#### Repeats

DifferentialAnalysis(data            = data_Repeats,
                     annotations     = annotations,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0_Repeats_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

DifferentialAnalysis(data            = data_Repeats,
                     annotations     = annotations_A1,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A1_Repeats_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

DifferentialAnalysis(data            = data_Repeats,
                     annotations     = annotations_A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2_Repeats_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

DifferentialAnalysis(data            = data_Repeats,
                     annotations     = annotations_AVUS,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUS_Repeats_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )

DifferentialAnalysis(data            = data_Repeats,
                     annotations     = annotations_AVUSA2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2_Repeats_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


DifferentialAnalysis(data            = data_Repeats,
                     annotations     = annotations_A0A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2_Repeats_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )




quit(save="no")
#stop("It only remains MixOmix and statistical codes")




anno <- annotations_A0A2[annotations_A0A2$Sample_Group == "ATM", ]
anno[anno$ClinVar %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic"), ]$ClinVar <- "PV"
anno[anno$ClinVar == "Uncertain significance", ]$ClinVar <- "VUS"


DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "ClinVar",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "PV - VUS",
                     comparison_text = "Comparison of PV and VUS ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2-PV_VUS_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )




# Comparison of Monash to monash samples

anno <- annotations_A0A2[annotations_A0A2$Source == "Monash", ]
#anno[anno$ClinVar %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic"), ]$ClinVar <- "PV"
#anno[anno$ClinVar == "Uncertain significance", ]$ClinVar <- "VUS"


DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non ATM cancers from MONASH",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2-Monash_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )



anno <- annotations_A2[annotations_A2$Source == "Monash", ]

DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non ATM cancers from MONASH",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2-Monash_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


anno <- annotations_AVUSA2[annotations_AVUSA2$Source == "Monash", ]

DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non ATM cancers from MONASH",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2-Monash_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )







# Comparison of only CoF-AT and GENESIS samples to Monash controls

anno <- annotations_A0A2[
    (annotations_A0A2$Sample_Group == "Non.ATM" & annotations_A0A2$Source == "Monash") |
    (annotations_A0A2$Sample_Group == "ATM" & annotations_A0A2$Source == "Curie"), ]

DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non ATM cancers from MONASH",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2-Curie_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


anno <- annotations_A2[
    (annotations_A2$Sample_Group == "Non.ATM" & annotations_A2$Source == "Monash") |
    (annotations_A2$Sample_Group == "ATM" & annotations_A2$Source == "Curie"), ]

DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non ATM cancers from MONASH",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2-Curie_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


anno <- annotations_AVUSA2[
    (annotations_AVUSA2$Sample_Group == "Non.ATM" & annotations_AVUSA2$Source == "Monash") |
    (annotations_AVUSA2$Sample_Group == "ATM" & annotations_AVUSA2$Source == "Curie"), ]

DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non ATM cancers from MONASH",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2-Curie_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )







# Test heterogeneity

# Comparison of ATM CoF-AT and GENESIS samples to ATM Monash samples

anno <- annotations_A0A2[
    (annotations_A0A2$Sample_Group == "ATM" & annotations_A0A2$Source == "Monash") |
    (annotations_A0A2$Sample_Group == "ATM" & annotations_A0A2$Source == "Curie"), ]

DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "Source",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "Curie - Monash",
                     comparison_text = "Comparison of ATM and Non ATM cancers from MONASH",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2-ATM-CurievsMonash_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


anno <- annotations_A2[
    (annotations_A2$Sample_Group == "ATM" & annotations_A2$Source == "Monash") |
    (annotations_A2$Sample_Group == "ATM" & annotations_A2$Source == "Curie"), ]

DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "Source",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "Curie - Monash",
                     comparison_text = "Comparison of ATM and Non ATM cancers from MONASH",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A2-ATM-CurievsMonash_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )


anno <- annotations_AVUSA2[
    (annotations_AVUSA2$Sample_Group == "ATM" & annotations_AVUSA2$Source == "Monash") |
    (annotations_AVUSA2$Sample_Group == "ATM" & annotations_AVUSA2$Source == "Curie"), ]

DifferentialAnalysis(data            = data_Promoters,
                     annotations     = anno,
                     colInterest     = "Source",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "Curie - Monash",
                     comparison_text = "Comparison of ATM and Non ATM cancers from MONASH",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2-ATM-CurievsMonash_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )










































library(ggpubr)

data = data_Promoters["CTLA4", anno$Sample_Name]
data <- t(data)
data <- merge(anno["ClinVar"], data, by="row.names")


comparisons <- list( c("PV", "VUS") )
t.test <- compare_means(CTLA4 ~ ClinVar, comparisons = comparisons,
                        p.adj= "BH", method='wilcox.test', data = data)

t.text <- t.test %>%
    mutate(y.position = c(-1.5))
t.test$label <- paste0(t.test$p.adj, " (", t.test$p.signif, ")")
print(t.test)


    # custom theme
    custom.theme <- theme(legend.position="none",
                          panel.grid = element_line(colour = "grey92"),
                          panel.grid.minor = element_line(linewidth = rel(0.5)),
                          panel.grid.major = element_line(linewidth = rel(0.5)))


data <- data %>% group_by(ClinVar) %>% mutate(count=n()) %>% data.frame()
        #data$Sample_Group <- paste0(data$Sample_Group, "\n(", data$count, ")")
        #print(data)
        
        p1 <- ggviolin(data,
                       x        = "ClinVar",
                       y        = "CTLA4",
                       #palette  = col, #c("#01BEC5", "#F8756E"),
                       fill     = "ClinVar",
                       add      = "boxplot") +
            #labs(title = paste0(main),
            #     subtitle = paste0("(", nprobes, " probes)")) +
            theme(legend.position="none",
                  plot.title = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size=rel(0.95))) + 
    #stat_compare_means(comparisons  = comparisons, method="t.test")# +
        stat_pvalue_manual(t.test, label = "label", y.position = 2) +
        xlab("Type of tumours") +
        custom.theme














# In Promoters
DifferentialAnalysis(data            = data_Promoters_Beta,
                     annotations     = annotations_A0A2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "A0A2-Beta_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )



DifferentialAnalysis(data            = data_Promoters_Beta,
                     annotations     = annotations_AVUSA2,
                     colInterest     = "Sample_Group",
                     analysisLevel   = "GenomeWide", #Gene or GenomeWide
                     #listProbes      = ProbesInATMGene ,
                     contrastToUse   = "ATM - Non.ATM",
                     comparison_text = "Comparison of ATM and Non.ATM cancers",
                     FC.threshold    = 1,
                     p.val.threshold = 0.05,
                     path            = file.path(OutFiles, "AVUSA2-Beta_Promoters_"),
                     GO              = FALSE,
                     KEGG            = FALSE,
                     keepAll         = TRUE
                     )



