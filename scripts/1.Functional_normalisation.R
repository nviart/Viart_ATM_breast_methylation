#!/usr/bin/env Rscript
#Torque Configuration
#PBS -l walltime=20:10:30
#PBS -l mem=100gb
#PBS -l nodes=2:ppn=4
#PBS -q batch
#PBS -N Normalisation_script
#PBS -j oe


#---- Script du 11.10.2021 ----#
#---- Script minimalisé pour normaliser les données EPIC avec la normalisation Functional (après benchmark des différentes méthodes de normalisation) ----#

# Nécessite le répertoire avec les données EPIC et le nom du fichier descriptif (dans le même dossier)
#EpicFolder = "/data/users/nviart/ATM_Analysis/data/methylation"
#DescriptorFileName1 = "Annotations_30.11.23_Curie.csv"
#DescriptorFileName2 = "Annotations_13.02.2024_Australia_ABCFR_MCCS.csv"

# Loading of customised functions and path to data to use
source(file.path(getwd(), "config.R"))

# Activate renv project
renv::activate(project = renv.path)


# Nécessite les fichiers spécifiant les sondes à supprimer
# Nécessite de spécification le chemin de sortie du fichier avec les valeurs Beta:
#library("optparse")
#option_list = list(
#    make_option(c("-d", "--outDir"), type="character", default=NULL,
#                help="output directory", metavar="character")
#)
#opt_parser = OptionParser(option_list=option_list)
#opt = parse_args(opt_parser)
#Origin = opt$outDir


# Define the directory where files will be saved
BetaOutDir = file.path(Origin, "matrices/")


if(!dir.exists(BetaOutDir))
{
    dir.create(BetaOutDir, recursive = TRUE)
}

# Chargement des fonctions customisées
#source(file="/data/users/nviart/ATM_Analysis/svn/analyse/script/methylation/functions.R")

# Chargement des packages
library("knitr")
library("rmdformats")
library("bumphunter")
library("RColorBrewer")
library("plyr")
library("dplyr")
library("Biobase")
library("data.table")
library("org.Hs.eg.db")
library("GenomeInfoDb")
library("RPMM")
library("lumi")
library("ENmix")
library("minfi")


# Modification du répertoire de travail (où les données EPIC sont présentes)
setwd(EpicFolder)


# Charger les données (descriptors et fichiers IDAT)
#targets <- read.metharray.sheet(EpicFolder, pattern = DescriptorFileName)

# Lire les données brutes à partir des fichiers IDAT
#rgSet <- read.metharray.exp(targets=targets, force=TRUE)

# give the samples descriptive names
#targets$ID <- paste(targets$Sample <- Name)
#sampleNames(rgSet) <- targets$ID
#rgSet

# Phenotypic data
#phenoData <- pData(rgSet)

# Manifest => probe design
#manifest <- getManifest(rgSet)



ReadData <- function(folder, DescriptorFileName,
                     differentVersion = FALSE,
                     whichVersion     = NULL)
{
    # Charger les données (descriptors et fichiers IDAT)
    targets <- read.metharray.sheet(EpicFolder, pattern = DescriptorFileName)
    targets$Basename <- gsub('"', "", targets$Basename)
    targets$Basename <- gsub('^c\\(|\\)$', "", targets$Basename)
    print(head(targets))

    if ("Selection" %in% colnames(targets)){
        targets <- targets[targets$Selection == "1", ]
    }
    
    if (differentVersion == TRUE){
        if(is.null(whichVersion)){
            stop("You indicated the use of different versions of arrays, please provide the type of array for which you want lo load the data")
        } else{
            targets <- targets[targets$ArrayVersion == whichVersion, ]
        }
    }
    #print(head(targets))
    # If somes files are not found, report them and suppress them from the targets
    if ("character(0" %in% targets$Basename){
        toPrint <- targets[targets$Basename == "character(0", ]
        print(toPrint)
        rm(toPrint)
        targets <- targets[targets$Basename != "character(0", ]
    }

    # For testing, select only the three first samples of each dataset
    #targets <- targets[1:4, ]
    #print(targets)
    
    # Lire les données brutes à partir des fichiers IDAT
    rgSet <- read.metharray.exp(targets=targets, force=TRUE)
    
    # give the samples descriptive names
    if ("Sample_Name" %in% names(targets)){
        targets$ID <- paste(targets$Sample_Name)
        sampleNames(rgSet) <- targets$ID
    }
    rgSet

    return(rgSet)
}


# Use of sink to write preprocessing results in a file (number of probes, number of samples)
sink(file.path(BetaOutDir, "stat_preprocessing.txt"))

if (TCGA == FALSE){
    if (onlyEPIC){
        # Read individually the different data
        ## ATM samples first
        data1 <- ReadData(EpicFolder, DescriptorFileName1)
  
        ## And australian samples (EPIC)
        data2.EPIC <- ReadData(file.path(EpicFolder, "Monash"), DescriptorFileName2,
                               differentVersion = TRUE, whichVersion = "EPIC")
    
        print("Dimentions of the original data (Nbr probes, Nbr samples)\n")
        cat("data1: ", dim(data1), "\n")
        cat("data2.EPIC: ", dim(data2.EPIC), "\n")
        
        # Combine all arrays
        rgSet <- combineArrays(object1=data1, object2=data2.EPIC)

    } else {
        # Read individually the different data
        ## ATM samples first
        data1 <- ReadData(EpicFolder, DescriptorFileName1)
        ## Then australian samples (450k)
        #data2.450k <- ReadData(file.path(EpicFolder, "Monash", "ABCFR_Breast_FFPE"), DescriptorFileName2,
        #                       differentVersion = TRUE, whichVersion = "HM450K")

        data2.450k <- ReadData(EpicFolder, DescriptorFileName2,
                               differentVersion = TRUE, whichVersion = "HM450K")
        
        ## And australian samples (EPIC)
        data2.EPIC <- ReadData(file.path(EpicFolder, "Monash"), DescriptorFileName2,
                               differentVersion = TRUE, whichVersion = "EPIC")

        print("Dimentions of the original data (Nbr probes, Nbr samples)\n")
        cat("data1: ", dim(data1), "\n")
        cat("data2.450k: ", dim(data2.450k), "\n")
        cat("data2.EPIC: ", dim(data2.EPIC), "\n")
        
        # Combine all arrays
        rgSet <- combineArrays(object1=data2.450k, object2=data2.EPIC)
        rgSet <- combineArrays(object1=rgSet, object2=data1)
    }
} else{ # SI TCGA data
    print("TCGA true")
    data = read.csv("/data/users/nviart/ATM_Analysis/data/GDCdata/TCGA-BRCA_variants.csv", header=T, sep=",")

    data$Patient_ID <- substr(data$barcode, 0, 12)

    # Get IDs to remove (variants in actionable genes)
    idsToRemove = unique(data[data$Selection == 0,]$Patient_ID)
    
    # Remove variants in actionable genes from data
    data = data[! data$Selection == 0,]

    # Get data to add ATM variant classification
    ATMclassif <- data[data$gene == "ATM", ]
    ATMclassif <- ATMclassif[c("Patient_ID", "ATM_classif")]
    # Remove duplicated rows
    ATMclassif <- unique(ATMclassif[, c("Patient_ID", "ATM_classif")])
    names(ATMclassif) <- c("Patient_ID", "ClinVar")

    table(ATMclassif$ClinVar)


    # Get identifiers of ATM variant carriers
    idsATM <- data[data$gene_id == "ATM",]$Patient_ID

    # 1218 fichiers téléchargés avec le manifest de la méthylation
    # Correspond bien aux 1218 fichiers présents dans le manifest (resFinal)
    resFinal <- read.csv("~/programmes/GDC/Query_manifest_methylation_constit_ERpositive.txt", sep="\t")
    names(resFinal)[names(resFinal) == "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"] <- "barcode"

    print(dim(resFinal))

    # Get location of methylation data
    resFinal$Basename <- paste0("/data/users/nviart/ATM_Analysis/data/GDCdata/raw_data/methylation/IDATfiles/", resFinal$filename)
    # Remove green and red information
    #resFinal$Basename <- replace(resFinal$Basename, "_Red.idat", "")
    resFinal$Basename <- gsub("_Red.idat", "", resFinal$Basename)
    resFinal$Basename <- gsub("_Grn.idat", "", resFinal$Basename)
    # Remove duplicated file names due to the removal of the color
    resFinal <- resFinal[! duplicated(resFinal$Basename),]

    dim(resFinal)
    # => Remains 609 samples

    # Remove Solid Tissue Normal and metastatic tumours
    #resFinal$Sample.Type <- substr(resFinal$barcode, 14, 15) 
    resFinal <- resFinal[! resFinal$Sample_Type %in% c("Metastatic", "Solid Tissue Normal"), ]

    dim(resFinal)
    # => Remains 539 samples

    # Remove participant with variant in other genes than ATM
    resFinal <- resFinal[! resFinal$bcr_patient_barcode %in% idsToRemove,]

    dim(resFinal)
    # => Remains 508 samples

    # Indicate if sample is from an ATM variant carrier or not
    resFinal$Sample_Group <- ifelse(resFinal$Patient_ID %in% idsATM, "ATM", "Non.ATM")

    table(resFinal$Sample_Group)
    # => 30 ATM and 478 sporadic samples
    
    print(dim(resFinal))
    
    # Add ATM classification to the table
    resFinal <- merge(resFinal, ATMclassif, by = "Patient_ID", all.x=T)
    
    table(resFinal$ClinVar)
    
    resFinal$ClinVar[resFinal$ClinVar == "Pathogenic/Likely_pathogenic"] <- "Pathogenic/Likely pathogenic"
    resFinal$ClinVar[is.na(resFinal$ClinVar)] <- ""
    
    # => 1 Likely pathogenic, 7 Pathogenic, 1 Pathogenic/Likely pathogenic and 21 VUS
    
    
    table(resFinal$Sample_Group, resFinal$ethnicity)
    
    table(resFinal$ClinVar, resFinal$ethnicity)
    
    table(resFinal$Sample_Group, resFinal$race)
    
    table(resFinal$ClinVar, resFinal$race)

    nATM <- nrow(resFinal[resFinal$Sample_Group == "ATM",])
    resFinal <- rbind(dplyr::sample_n(resFinal[resFinal$Sample_Group == "Non.ATM",], nATM),
                      resFinal[resFinal$Sample_Group == "ATM",])
    resFinal$Study <- "TCGA-BRCA"
    
    # Read data
    rgSet <- read.metharray.exp("/data/users/nviart/ATM_Analysis/data/GDCdata/raw_data/methylation/IDATfiles", recursive=TRUE, targets=resFinal)

    # Modify names of the rgSet
    annotations <- data.frame(rgSet@colData)
    annotations$noid <- colnames(rgSet)
    colnames(rgSet) <- annotations[colnames(rgSet), ]$barcode
}


print("Dimentions of the combined arrays (Nbr probes, Nbr samples)\n")
cat("data combined: ", dim(rgSet), "\n")

print("table [Sample_Group, Study]: ")
print(table(pData(rgSet)$Sample_Group, pData(rgSet)$Study))

# Phenotypic data
phenoData <- pData(rgSet)

# Manifest => probe design
manifest <- getManifest(rgSet)

print("manifest: ")
print(manifest)

#system("wget https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
# system("unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")

library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)

# Force Minfi to use the new package
#rgSet@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
#print("forced")
#annotation(rgSet)["annotation"] = "ilm10b5.hg38"
#annotation(rgSet)["array"] = "IlluminaHumanMethylationEPIC"

# Load the new manifest file
mf=manifest.path
annEPIC <- readmanifest(mf)
annEPIC <- annEPIC$assay

write.table(annEPIC, file = file.path(BetaOutDir, 'manifest.tsv'), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

print("manifest written")

# Write the localisation of the probes in a file to map them with genes and promoters later
df_probes_locations <- annEPIC[, c("CHR_hg38", "Start_hg38", "End_hg38", "Name", "Strand_hg38")]
df_probes_locations$Start_hg38 <- as.numeric(df_probes_locations$Start_hg38)
df_probes_locations <- df_probes_locations[order(df_probes_locations$CHR_hg38, df_probes_locations$Start_hg38), ]
df_probes_locations <- df_probes_locations[complete.cases(df_probes_locations$CHR_hg38), ]
df_probes_locations <- df_probes_locations[!duplicated(df_probes_locations), ]
write.table(df_probes_locations, file = '/data/users/nviart/ATM_Analysis/data/annotation/updateFinal/Probes_annotations.bed', col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

print("positions writen")

#anno <- annEPIC$manifest
#manifestList <- annEPIC$manifestList

#maniTmp <- minfi:::read.manifest.EPIC(mf)
#anno <- maniTmp$manifest
#manifestList <- maniTmp$manifestList

# Manifest package
#IlluminaHumanMethylationEPICmanifest <- do.call(IlluminaMethylationManifest,
#                                                list(TypeI = manifestList$TypeI,
#                                                     TypeII = manifestList$TypeII,
#                                                     TypeControl = manifestList$TypeControl,
#                                                     TypeSnpI = manifestList$TypeSnpI,
#                                                     TypeSnpII = manifestList$TypeSnpII,
#                                                     annotation = "IlluminaHumanMethylationEPIC"))

#usethis::use_data(IlluminaHumanMethylationEPICmanifest, internal=TRUE)

#rgSet@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")

#rm(anno, manifestList)

#saveRDS(rgSet, file=file.path(BetaOutDir, "RGSET_BeforeFirstFiltering.rds"))
#rgSet <- readRDS(file.path(BetaOutDir, "RGSET_BeforeFirstFiltering.rds"))


# Filtrage des données 
## Etape 1
#detP <- minfi::detectionP(rgSet) > 0.01
detP <- minfi::detectionP(rgSet) > 0.01
keep <- colMeans(detP) < 0.10
rgSet <- rgSet[,keep]

#saveRDS(rgSet, file=file.path(BetaOutDir, "RGSET_FirstFiltering2.rds"))
phenoData <- subset(phenoData, rownames(phenoData) %in% colnames(rgSet))

print("Dimentions of arrays after first filtering(Nbr probes, Nbr samples)\n")
cat("data combined: ", dim(rgSet), "\n")

print("table [Sample_Group, Study]: ")
print(table(pData(rgSet)$Sample_Group, pData(rgSet)$Study))


## Etape 2
#hp = detectionP(rgSet) > 0.01
hp = detectionP(rgSet[, pData(rgSet)[pData(rgSet)$Sample_Group == "Sporadic", ]$Sample_Name]) > 0.01
exclude.hpv = rownames(hp)[rowMeans(hp) > 0.10]

rgSet = subsetByLoci(
            rgSet = rgSet,
            excludeLoci = c(exclude.hpv))
rm(exclude.hpv, keep, detP, hp)

print("Dimentions of arrays after second filtering(Nbr probes, Nbr samples)\n")
cat("data combined: ", dim(rgSet), "\n")

print("table [Sample_Group, Study]: ")
print(table(pData(rgSet)$Sample_Group, pData(rgSet)$Study))


if (probeFiltering)
{
    # Filtrage des sondes avec des SNP et dans des régions cross-réactives
    ## cross-reactive/non-specific
    cross.react <- read.csv(file.path(path.filters, '48639-non-specific-probes-Illumina450k.csv'), head = T, as.is = T)
    
    cross.react.probes <- as.character(cross.react$TargetID)

    ## BOWTIE2 multi-mapped
    multi.map <- read.csv(file.path(path.filter, 'HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt'), head = F, as.is = T)
    multi.map.probes <- as.character(multi.map$V1)

    ## determine unique probes
    filter.probes <- unique(c(cross.react.probes, multi.map.probes))

    ## filter rgSet
    rgSet = subsetByLoci(
        rgSet = rgSet,
        excludeLoci = c(filter.probes))

    print("Dimentions of arrays after removal of cross-reactive probes and multi-mapped probes (Nbr probes, Nbr samples)\n")
    cat("data combined: ", dim(rgSet), "\n")

    print("table [Sample_Group, Study]: ")
    print(table(pData(rgSet)$Sample_Group, pData(rgSet)$Study))

    ## remove data from memory
    rm(multi.map, multi.map.probes, cross.react, cross.react.probes, filter.probes)

    ## probes from Pidsley 2016 (EPIC)
    epic.cross1 <- read.csv(file.path(path.filters, '13059_2016_1066_MOESM1_ESM.csv'), head = T)
    epic.variants1 <- read.csv(file.path(path.filters, '13059_2016_1066_MOESM4_ESM.csv'), head = T)
    epic.variants2 <- read.csv(file.path(path.filters, '13059_2016_1066_MOESM5_ESM.csv'), head = T)
    epic.variants3 <- read.csv(file.path(path.filters, '13059_2016_1066_MOESM6_ESM.csv'), head = T)
    ## additional filter probes
    epic.add.probes <- c(as.character(epic.cross1$X), as.character(epic.variants1$PROBE), as.character(epic.variants2$PROBE),
                         as.character(epic.variants3$PROBE))
    ## final list of unique probes
    epic.add.probes <- unique(epic.add.probes)

    ## filter rgSet
    rgSet = subsetByLoci(
        rgSet = rgSet,
        excludeLoci = c(epic.add.probes))

    print("Dimentions of arrays after removal of probes from Pidsley et al (Nbr probes, Nbr samples)")
    cat("data combined: ", dim(rgSet), "\n")

    print("table [Sample_Group, Study]: ")
    print(table(pData(rgSet)$Sample_Group, pData(rgSet)$Study))


    ## remove data from memory
    rm(epic.cross1, epic.variants1, epic.variants2, epic.variants3, epic.add.probes)
}

# Save the object of read data
save(rgSet, file = file.path(BetaOutDir, "rgSet.RData"))


# system("wget https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
# system("unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
#mf="/data/users/nviart/ATM_Analysis/data/methylation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"
#annEPIC <- readmanifest(mf)
#annEPIC <- annEPIC$assay

#keep <- !(featureNames(rgSet) %in% annEPIC$Name[annEPIC$CHR_hg38 %in% c("chrX","chrY")])
#print(length(keep))
#rgSet <- rgSet[keep,]


#sex.probes <- annEPIC$Name[annEPIC$CHR_hg38 %in% c("chrX","chrY")]
#rgSet = subsetByLoci(
#            rgSet = rgSet,
#            excludeLoci = c(sex.probes))
#print(sex.probes)


# Installation of the IlluminaHumanMethylationEPICanno.ilm10b5.hg38
#(installed the 11.03.2022)
#install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest")
#install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")


## Normalisation des données
MSet_Funnorm <- preprocessFunnorm(rgSet)

print("Dimentions of arrays after normalisation (Nbr probes, Nbr samples)\n")
cat("data combined: ", dim(rgSet), "\n")

print("table [Sample_Group, Study]: ")
print(table(pData(rgSet)$Sample_Group, pData(rgSet)$Study))


# Sélection des sondes dans les chromosomes sexuels
sex.probes <- annEPIC$Name[annEPIC$CHR_hg38 %in% c("chrX","chrY")]

# Suppression
MSet_Funnorm <- MSet_Funnorm[!rownames(MSet_Funnorm) %in% sex.probes, ]

print("Dimentions of arrays after removal probes mapped to sexual chromosomes (Nbr probes, Nbr samples)\n")
cat("data combined: ", dim(rgSet), "\n")

print("table [Sample_Group, Study]: ")
print(table(pData(rgSet)$Sample_Group, pData(rgSet)$Study))

# Mappage sur le génome
gMSet_Funnorm <- mapToGenome(MSet_Funnorm, mergeManifest = TRUE)
rm(MSet_Funnorm)



## Suppress samples which are not anymore required
# Suppress the duplicated sample 
#sample.names = c('T0015dup')
#gMSet_Funnorm <- gMSet_Funnorm[, !colnames(gMSet_Funnorm) %in% sample.names]

# Suppress the non-cancerous samples
#sample.names = subset(phenoData, phenoData$Sample_Group == "Non-cancerous")$Sample_Name
#gMSet_Funnorm <- gMSet_Funnorm[, !colnames(gMSet_Funnorm) %in% sample.names]

# Subset annotations to have only remaining samples
gMSet_Funnorm@colData <- gMSet_Funnorm@colData[colnames(gMSet_Funnorm),]

# Save the raw normalized data
save(gMSet_Funnorm, file = file.path(BetaOutDir, "MethylSet_Funnorm_Normalisation.RData"))


# Sauvegarde des données normalisées
all_beta <- data.frame(getBeta(gMSet_Funnorm))
write.csv(all_beta, file.path(BetaOutDir, "Beta_Funnorm_Normalisation.csv"))
save(all_beta, file = file.path(BetaOutDir, "Beta_Funnorm_Normalisation.RData"))

print("Dimentions of Beta matrix (Nbr probes, Nbr samples)\n")
cat("data combined: ", dim(all_beta), "\n")

getMCustom <- function(Beta){
       M <- log2(Beta / (1 - Beta))
       return(M)
}


# Calcul des M values
all_M_Funnorm <- getMCustom(all_beta)
# Filtre des valeurs infinies ou manquantes
all_M_Funnorm[all_M_Funnorm == "Inf"] <- NA
all_M_Funnorm[all_M_Funnorm == "-Inf"] <- NA
all_M_Funnorm <- all_M_Funnorm[complete.cases(all_M_Funnorm), ]

#all_M_Funnorm <- data.frame(getM(gMSet_Funnorm))
all_M <- data.frame(all_M_Funnorm)
write.csv(all_M, file.path(BetaOutDir, "M_Funnorm_Normalisation.csv"))
save(all_M, file = file.path(BetaOutDir, "M_Funnorm_Normalisation.RData"))

print("Dimentions of M matrix (Nbr probes, Nbr samples)\n")
cat("data combined: ", dim(all_M), "\n")

annotations <- data.frame(gMSet_Funnorm@colData)
annotations$Sample_Group[annotations$Sample_Group == "Sporadic"] <- "Non.ATM"
annotations$ClinVar[annotations$ClinVar == "Pathogenic ; Pathogenic"] <- "Pathogenic"
annotations$ClinVar[annotations$ClinVar == "Pathogenic*"] <- "Pathogenic"


if (TCGA == FALSE){
    ## !!!!!! Parce qu'il n'y a que des échantillons AT dans CoFAT pour le moment
    HBOCFamilies <- annotations$Sample_Name[annotations$Study == "GENESIS"]
    ATFamilies <- annotations$Sample_Name[annotations$Study == "CoFAT"]
    annotations$FamilySelection <- ifelse(annotations$Sample_Name %in% HBOCFamilies, "HBOC",
                                          ifelse(annotations$Sample_Name %in% ATFamilies, "AT", "Non variant carriers")
                                          )
}

write.csv(annotations, file.path(BetaOutDir, "Annotations.csv"))

# Pour charger le fichier, utiliser la commande "read.csv(file.path(BetaOutDir, "Beta_Funnorm_Normalisation.csv"), row.names = 1)"
