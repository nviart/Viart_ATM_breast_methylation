#!/usr/bin/env Rscript
#Torque Configuration
#PBS -l walltime=20:10:30
#PBS -l mem=60gb
#PBS -l nodes=2:ppn=4
#PBS -q batch
#PBS -N Exploratory_analysis
#PBS -j oe

source("~/ATM_Analysis/svn/analyse/script/methylation/MethPipeline/config.R")

# Activate renv project
renv::activate(project = renv.path)

# Lien vers les annotations initiales (avec tous les échantillons)
InitialAnnotationFileInitialAnnotationFile = "~/ATM_Analysis/data/methylation/Annotations_30.11.23_Curie.csv"
# Lien vers les annotations et les Beta values
#library("optparse")
#option_list = list(
#    make_option(c("-d", "--outDir"), type="character", default=NULL,
#                help="output directory", metavar="character")
#)
#opt_parser = OptionParser(option_list=option_list)
#opt = parse_args(opt_parser)

#Origin = file.path(opt$outDir)

#Origin = "/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation/"
DirFiles = file.path(Origin, "matrices/")
# Lien vers le dossier de sortie
OutFilesF = file.path(Origin, "paper/")

if(!dir.exists(OutFilesF))
{
    dir.create(OutFilesF, recursive = TRUE)
}

# Chargement des librairies
library("RColorBrewer")
library("minfi")
library('dplyr')
library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
library("ggplot2")
library("org.Hs.eg.db")
library("GenomeInfoDb")
library("cowplot")
library("scales")
library("dendextend")
library("ENmix")
library("ggpubr")
library("readxl")
library("tidyr")

# Lecture du fichier des Beta values
#all_beta <- read.csv(file.path(DirFiles, "Beta_Funnorm_Normalisation.csv"), row.names = 1)
#load(file.path(DirFiles, "Beta_Funnorm_Normalisation.RData"))

# Lecture du fichier des M values (Les M values ayant été calculées dans le script précédent et les valeurs infinies supprimées)
#all_M <- read.csv(file.path(DirFiles, "M_Funnorm_Normalisation.csv"), row.names = 1)
load(file.path(DirFiles, "M_Funnorm_Normalisation.RData"))

# Chrgement des données mappées
load(file.path(DirFiles, "Mapping.RData"))

# chargement des annotations
annotations <- read.csv(file.path(DirFiles, "Annotations.csv"), row.names = 1)

# Add LOH information
ATM.LOH.Yes <- c("T0220", "T0076", "T0099", "T0123", "T0248", "T0077R", "T0072", "T0078")
ATM.LOH.No <- c("T0120", "T0015", "T0009", "T0249", "T0247")    
annotations[["ATM_LOH"]]  <- ifelse(annotations$Sample_Name %in% ATM.LOH.Yes, "LOH",
                             ifelse(annotations$Sample_Name %in% ATM.LOH.No, "No LOH",
                                    "Unknown LOH")
                             )

# Modify variant
annotations[annotations$Sample_Name == "S2013401", ]$ClinVar <- "Pathogenic"
# Modify Date diag
annotations$Date_diag <- stringr::str_sub(annotations$tum_date_diag, -4, -1)
annotations$Date_diag <- as.numeric(annotations$Date_diag)
# Modiy ER status
annotations$ER_status <- as.character(annotations$ER_status)

# Chargement des fonctions customisées
source(file = file.path(function.path, "functions.R"))

print("data loaded")

# Exploratory Analysis
# For this exploratory analysis, the 10000 probes having the highest inter-quartile range were selected.
# Sélection du nombre désiré de sondes les plus variables
n_most_variable <- 10000

# Sélection des sondes les plus variables en utilisant les Beta values
#pIQBeta <- variable_feature_selection(all_beta, thres.num = n_most_variable)

# Sélection des sondes les plus variables en utilisant les M values
pIQM <- variable_feature_selection(all_M, thres.num = n_most_variable)





##############################################################
###################### Figure S3 and S4 #######################
###############################################################

annotations$ER_status <- as.character(annotations$ER_status)
#annotations[annotations$Variant_type == "splice acceptor",]$Variant_type_aggregated <- "LoF"
annotations$ER_status2 <- annotations$ER_status
annotations[is.na(annotations$ER_status2), ]$ER_status2 <- "NA"
annotations[annotations$ER_status2 == "1", ]$ER_status2 <- "Positive"
annotations[annotations$ER_status2 == "0", ]$ER_status2 <- "Negative"
annotations[annotations$ER_status2 == "NA", ]$ER_status2 <- NA

annotations$ClinVar2 <- ifelse(annotations$ClinVar %in% c("Likely pathogenic",
                                                          "Pathogenic",
                                                          "Pathogenic/Likely pathogenic"), "PV",
                        ifelse(annotations$ClinVar %in% c("Uncertain significance"), "VUS", "Non-carrier"))

#umap <- umap_func_publi(matrix                      = all_M[pIQM, ],
#                        sampleTable                 = annotations,
#                        highlightedVar              = "ER_status2",
#                        highlightedVarName          = "ER status",
#                        sampleLabel                 = "Sample_Name",
#                        sampleLabelColumnCondition  = "Sample_Group",
#                        sampleLabelValueCondition   = "ATM",
#                        shape                       = "ClinVar2",
#                        colors                      = brewer.pal(n = 3, name = "Dark2"),
#                        fontsize                    = 14,
#                        pointSize                   = 2,
#                        metric                      = "pearson2",
#                        title                       = "UMAP using the 10.000 most variable probes") 

#ggsave(file.path(OutFiles, "UMAP_publication.png"), umap, height=7, width=12)


annotations$Group <- ifelse(annotations$Sample_Group == "ATM", "ATM variant carriers",
                            "Non-carriers")

umap <- umap_func_publiv2(matrix                      = all_M[pIQM, ],
                          sampleTable                 = annotations,
                          highlightedVar              = "ER_status2",
                          highlightedVarName          = "ER status",
                          #sampleLabel                 = "Sample_Name",
#                        sampleLabelColumnCondition  = "Sample_Group",
#                        sampleLabelValueCondition   = "ATM",
                          shape                       = "Group",
                          colors                      = brewer.pal(n = 3, name = "Dark2"),
                          fontsize                    = 14,
                          pointSize                   = 2,
                          metric                      = "pearson2",
                          title                       = NULL,#"UMAP using the 10.000 most variable probes",
                          colCond = "Sample_Group",
                          valCond = "ATM",
                          seed = 1234) 


ggsave(file.path(OutFilesF, "Figure S3.png"), umap, height=7, width=12, dpi=400)




annotations$Series <- ifelse(annotations$Source == "Curie", "French series",
                            ifelse(annotations$Source == "Monash", "Australian series", ""))
annotations$Series <- factor(annotations$Series, c("French series", "Australian series"))

umap <- umap_func_publiv2(matrix                      = all_M[pIQM, ],
                          sampleTable                 = annotations,
                          highlightedVar              = "Series",
                          highlightedVarName          = "Series",
                          shape                       = "Group",
                          colors                      = brewer.pal(n = 3, name = "Set1"),
                          fontsize                    = 14,
                          pointSize                   = 2,
                          metric                      = "pearson2",
                          title                       = NULL,#"UMAP using the 10.000 most variable probes",
                          colCond = "Sample_Group",
                          valCond = "ATM",
                          seed = 1234) 

ggsave(file.path(OutFilesF, "Figure S4.png"), umap, height=7, width=12, dpi=400)







### Statistical analysses
library(dplyr)

sink(file = file.path(OutFilesF, "statistics.txt"))

#annotations <- read.csv(file.path(DirFiles, "Annotations.csv"), row.names = 1)

#source(file="/data/users/nviart/ATM_Analysis/svn/analyse/script/methylation/functions.R")


# Modify groups
#annotations$group <- paste0(annotations$Study, "_", annotations$Sample_Group)
#annotations$group <- factor(annotations$group,
#                            levels = c("CoFAT_ATM", "GENESIS_ATM", "ABCFR_ATM", "ABCFR_Non.ATM", "MCCS_Non.ATM"))

# Modify ClinVar classification
annotations$ClinVar[annotations$ClinVar %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")] <- "PV"
annotations$ClinVar[annotations$ClinVar %in% c("Uncertain significance")] <- "VUS"


# Number of samples
table <- annotations %>% group_by(Sample_Group, ClinVar) %>%
    summarize(count = n())
table

# Diagnosis age
table <- annotations %>% group_by(Sample_Group, ClinVar) %>%
    summarize_at(vars(Diagnosis_age),
                 list(mean = mean, min = min, max = max, std=sd))
table

table <- annotations %>% group_by(Sample_Group) %>%
    summarize_at(vars(Diagnosis_age),
                 list(mean = mean, min = min, max = max, std=sd))
table

# Diagnosis year
table <- annotations %>% group_by(Sample_Group, ClinVar) %>%
    summarize_at(vars(Date_diag),
                 list(mean = mean, min = min, max = max, std=sd, med = median))
table %>% print(n=40, width=Inf)


table <- annotations %>% group_by(Sample_Group) %>%
    summarize_at(vars(Date_diag),
                 list(mean = mean, min = min, max = max, std=sd, med = median))
table %>% print(n=40, width=Inf)




# ER_status
annotations %>% group_by(Sample_Group, ClinVar) %>% count(ER_status)






#############
# ER status #
#############

# ER status
print(prop.test(x = c(23, 350), n = c(23+3, 350+107), alternative="greater"))

#########
# Grade #
#########

# Grade
print(prop.test(x = c(12+9, 171+186), n = c(12+9+3, 171+186+99), alternative="greater"))




sink()











#############################################################
######################### Figure S2 #########################
#############################################################

# Boxplot with only PV carriers

annotations_A1$group <- annotations_A1$Sample_Group
annotations_A1$group <- ifelse(annotations_A1$group == "Non.ATM", "Non-carriers", "ATM PV carriers") 
annotations_A1$group <- factor(annotations_A1$group, levels = c("ATM PV carriers", "Non-carriers"))


# custom theme

theme_Publication <- function(base_size=12) {#, base_family="Times") {
    library(grid)
    (
        theme_bw(base_size=base_size)#, base_family=base_family)
        + theme(#legend.position = "left",
              plot.title = element_text(hjust = 0.5),
              axis.text = element_text(size = base_size),
              panel.grid = element_line(colour = "grey92"),
              panel.grid.minor = element_line(linewidth = rel(0.5)),
              panel.grid.major = element_line(linewidth = rel(0.5))
          )
    )
}

     
require(rstatix)
require(ggpubr)


####################
# Age at diagnosis #
####################

print("Age at diagnosis")
print(tapply(annotations_A1$Diagnosis_age, annotations_A1$group, summary))

# Normality test
annotations_A1 %>% group_by(group) %>% shapiro_test(Diagnosis_age)

# Year of diagnosis are not normally distributed => Wilcoxon test

# test des variances
#var.test(Diagnosis_age ~ Sample_Group, data = annotations_A1)
# => p-value of 0.76 => no significant diffrenece between the two variances
comparisons <- list( c("Control", "ATM") )
t.test <- ggpubr::compare_means(Diagnosis_age ~ group, comparisons = comparisons,
                                p.adj= "BH", method='wilcox.test', data = annotations_A1,
                                paired=FALSE)
t.test <- t.test %>%
    mutate(y.position = c(-1.5))
t.test$label <- paste0("p=", t.test$p.adj)#, " (", t.test$p.signif, ")")
print(t.test)


plot1 <- ggpubr::ggviolin(annotations_A1,
                          x="group",
                          y="Diagnosis_age",
                          fill="group",
                          add = "boxplot",
                          palette = c("#B72B22", "#01BEC5")) +
theme_Publication() + 
ggtitle(paste0("Age at diagnosis")) +
    stat_pvalue_manual(t.test, label = "label", y.position = 85) +
    ylab("Age at diagnosis") +
    theme_Publication() + 
    theme(legend.position="none",
          plot.title = element_blank(),
          axis.title.x=element_blank()) +
    scale_x_discrete(labels=c("ATM heterozygous PV carriers", "Noncarriers"))

#ggsave(file.path(OutFilesF, "ViolinPlot_age_diagnosis_control.png"), plot1, dpi=300)


# UMAP visualisation
umap <- umap_func_publi(matrix                      = all_M[pIQM, ],
                        returnData                  = TRUE,
                        sampleTable                 = annotations,
                        highlightedVar              = "Diagnosis_age",
                        highlightedVarName          = "Age at diagnosis",
                        sampleLabel                 = "Sample_Name",
                        sampleLabelColumnCondition  = "group",
                        sampleLabelValueCondition   = "ATM",
                        fontsize                    = 14,
                        pointSize                   = 2,
                        metric                      = "pearson2",
                        title                       = "UMAP using the 10.000 most variable probes") 

g1 <- ggplot(data = umap,
             mapping = aes(x = UMAP1, y = UMAP2, color = Diagnosis_age)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
#    ggtitle("UMAP using the 10.000 most variable probes") +
    labs(color = "Age") + 
    theme_bw() +
    theme_Publication() + 
    theme(plot.title = element_blank()) + geom_point() +
    scale_x_discrete(labels=c("ATM heterozygous PV carriers", "Noncarriers"))






#####################
# Year of diagnosis #
#####################

print("Year of diagnosis")
print(tapply(annotations_A1$Date_diag, annotations_A1$group, summary))

# Year of diagnosis are not normally distributed => Wilcoxon test
# test des variances
#var.test(Date_diag ~ Sample_Group, data = annotations_A1)
# => p-value of 0.19 => no significant diffrenece between the two variances
comparisons <- list( c("Control", "ATM") )
t.test <- compare_means(Date_diag ~ group, comparisons = comparisons,
                        p.adj= "BH", method='wilcox.test', data = annotations_A1)
t.test <- t.test %>%
    mutate(y.position = c(-1.5))
t.test$label <- paste0("p=", t.test$p.adj)#, " (", t.test$p.signif, ")")
print(t.test)

plot2 <- ggviolin(annotations_A1,
                  x="group",
                  y="Date_diag",
                  fill="group",
                  add = "boxplot",
                  palette = c("#B72B22", "#01BEC5")) +
ggtitle(paste0("Year of diagnosis")) +
    stat_pvalue_manual(t.test, label = "label", y.position = 2020) +
    ylab("Year of diagnosis") + 
    theme_Publication() + 
    theme(legend.position="none",
          plot.title = element_blank(),
          axis.title.x=element_blank()) +
    scale_x_discrete(labels=c("ATM heterozygous PV carriers", "Noncarriers"))

#ggsave(file.path(OutFilesF, "ViolinPlot_date_diag_control.png"), plot2, dpi=300)



# UMAP visualisation
umap <- umap_func_publi(matrix                      = all_M[pIQM, ],
                        returnData                  = TRUE,
                        sampleTable                 = annotations,
                        highlightedVar              = "Date_diag",
                        highlightedVarName          = "Year of diagnosis",
                        sampleLabel                 = "Sample_Name",
                        sampleLabelColumnCondition  = "group",
                        sampleLabelValueCondition   = "ATM",
                        fontsize                    = 14,
                        pointSize                   = 2,
                        metric                      = "pearson2",
                        title                       = "UMAP using the 10.000 most variable probes") 

g2 <- ggplot(data = umap,
             mapping = aes(x = UMAP1, y = UMAP2, color = Date_diag)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") + 
    ggtitle("UMAP using the 10.000 most variable probes") +
    labs(color = "Year") + 
    theme_bw()+
    theme_Publication() + 
    theme(plot.title = element_blank()) + geom_point() +
    scale_x_discrete(labels=c("ATM heterozygous PV carriers", "Noncarriers"))



allPlots <- ggpubr::ggarrange(plot1, g1, plot2, g2, labels = c("A", "B", "C", "D"), ncol = 2, nrow=2)

ggsave(file.path(OutFilesF, "Figure S2.png"), allPlots,
       dpi=300, units="cm", width=45, height=35, scale = 0.7)













###################################################################
########################## Figure 1 ###############################
###################################################################


plot <- function(data             = M.Promoters.annotated,
                 Annotations      = annotations_A0A2,
                 Gene             = "ATM",
                 accordingTo      = NULL,
                 significance     = TRUE,
                 title            = NULL,
                 yposition        = NULL,
                 compareControls  = NULL,
                 xlab,
                 y1 = NULL,
                 y2 = NULL,
                 returnData = FALSE)
{
    require(ggpubr)
    require(dplyr)
    
    if (! "ATM_LOH" %in% colnames(Annotations))
    {
        stop("Provide ATM_LOH annotations")
    }
    # Create dataframe
    a <- data.frame(row.names    = rownames(Annotations),
                    Sample_Group = Annotations[["Sample_Group"]],
                    ClinVar      = Annotations[["ClinVar"]],
                    ATM_LOH      = Annotations[["ATM_LOH"]]#,
                    )

    a[[accordingTo]] <- Annotations[[accordingTo]]
    
    a[a$ClinVar %in% c("Pathogenic", "Likely pathogenic",
                       "Pathogenic/Likely pathogenic"),]$ClinVar <- "ATM PV"
    if ("Uncertain significance" %in% a$ClinVar)
    {
        a[a$ClinVar %in% c("Uncertain significance"),]$ClinVar <- "ATM VUS"
    }
    if ("" %in% a$ClinVar | NA %in% a$ClinVar)
    {
        a[a$ClinVar %in% c("", NA),]$ClinVar <- "Control"
    }
    a$ClinVar <- factor(a$ClinVar, levels = c("ATM PV", "ATM VUS", "Control"))
    
    print(head(a))
    a <- merge(a, t(data[Gene, ]), by = "row.names")

    if ("Inactive" %in% colnames(a)){
        a$Inactive <- factor(a$Inactive, levels = c("active ATM", "inactive ATM", "Control"))
    }

    print(unique(a$Sample_Group))
    a[a$Sample_Group == "Non.ATM", ]$Sample_Group <- "Control"
    a$Sample_Group <- factor(a$Sample_Group, levels = c("ATM", "Control"))
    
    a$ATM_LOH <- ifelse(a$Sample_Group == "Control", "Control",
                        paste0(a$Sample_Group, " ", a$ATM_LOH))
    
    a$ATM_LOH <- factor(a$ATM_LOH,
                        levels = c("ATM LOH", "ATM No LOH", "ATM Unknown LOH", "Control"))

    print(head(a))

    NameList <- c("Sample", "Sample_Group", "ClinVar", "ATM_LOH")
    if (accordingTo %in% c("Sample", "Sample_Group", "ClinVar", "ATM_LOH")){
        names(a) <- c(NameList,  "mean")
    } else {names(a) <- c(NameList, accordingTo, "mean")}
    

#    if (!accordingTo == "ATM_LOH"){
#        names(a) <- c("Sample", "Sample_Group", "ClinVar", "ATM_LOH", accordingTo, "mean")
#        } else {names(a) <- c("Sample", "Sample_Group", "ClinVar", accordingTo, "mean")}

    print(head(a))
    if (!is.null(compareControls))
    {
        controls <- a %>% filter(ClinVar == "Control") %>% dplyr::select(mean)
        meanControls <- mean(controls[,1])
        sdControls <- sd(controls[,1])
        print("mean for controls")
        print(meanControls)
        print("sd for controls")
        print(sdControls)
    }

    if (length(unique(a[[accordingTo]])) == 3){
        pal <- c("#B72B22", "#FF5349", "#01BEC5")
    } else {
        if (length(unique(a[[accordingTo]])) == 2){
            pal <- c("#B72B22", "#01BEC5")
        } else {
            pal <- c("#B72B22", "#FF5349", "#F8756D", "#01BEC5")
        }
    }
    
    if (accordingTo == "LOH")
    {
        xlab = "ATM LOH"
        if (significance)
        {
            max = max(a$mean) + 0.5
            if (is.null(yposition))
            {
                yposition <- c(max + 0.5 , max, max - 0.5)
            }
            
            comparisons <- list(c("ATM LOH", "Control"),
                                c("ATM No LOH", "Control"),
                                c("ATM Unknown LOH", "Control")
                                )
            t.test <- compare_means(mean ~ ATM_LOH, comparisons = comparisons,
                                    p.adj= "BH", method='wilcox.test', data = a)
            t.text <- t.test %>%
                mutate(y.position = max)
            t.test$label <- paste0(t.test$p.adj, " (", t.test$p.signif, ")")
            
            print(t.test)
            #t.test <- t.test %>% filter(group1 == "Control" | group2 == "Control")
        }
        if (returnData)
        {
            return(a)
        }
        p1 <- ggviolin(a,
                       x        = "ATM_LOH",
                       y        = "mean",
                       palette  = pal,
                       fill     = "ATM_LOH",
                       add      = "boxplot",
                       add.params = list(outlier.shape=0),
                       outlier.colour  = NA) +
            xlab(xlab)
    } else{
        print("step2")
        
        if (significance & (accordingTo == "ClinVar" | accordingTo == "ATM_LOH") ){
            if (accordingTo == "ClinVar"){
                max = max(a$mean) + 0.5
                if (is.null(yposition)){
                    yposition <- c(max + 0.5 , max, max - 0.5)
                }
                comparisons <- list(c("ATM PV ER+", "Control"),
                                    c("ATM VUS ER+", "Control"),
                                    c("ATM PV ER+", "ATM VUS ER+")
                                    )
                t.test <- compare_means(mean ~ ClinVar, comparisons = comparisons,
                                        p.adj= "BH", method='wilcox.test', data = a)
 
            }
            if (accordingTo == "ATM_LOH"){
                max = max(a$mean) + 0.5
                if (is.null(yposition)){
                    yposition <- c(max + 0.5 , max, max - 0.5)
                }
                comparisons <- list(c("ATM No LOH", "Control"),
                                    c("ATM Unknown LOH", "Control")
                                    )
                t.test <- compare_means(mean ~ ATM_LOH, comparisons = comparisons,
                                        p.adj= "BH", method='wilcox.test', data = a)
            }
            if (accordingTo == "Sample_Group"){
                max = max(a$mean) + 0.5
                if (is.null(yposition)){
                    yposition <- c(max + 0.5 , max, max - 0.5)
                }
                comparisons <- list(c("ATM", "Control"))
                t.test <- compare_means(mean ~ Sample_Group, comparisons = comparisons,
                                        p.adj= "BH", method='wilcox.test', data = a)   
            }
    
            t.text <- t.test %>%
                mutate(y.position = max)
            t.test$label <- paste0(t.test$p.adj, " (", t.test$p.signif, ")")
            a <- a %>% group_by(ClinVar) %>% mutate(count=n()) %>% data.frame()       
        } else {
            max = max(a$mean) + 0.5
            if (is.null(yposition)){
                yposition <- c(max + 0.5 , max, max - 0.5)
            }
            form <- paste0("mean", " ~ ", accordingTo)
            t.test <- compare_means(formula(form),
                                    p.adj= "BH", method='wilcox.test', data = a)
            t.text <- t.test %>%
                mutate(y.position = max)
            t.test$label <- paste0(t.test$p.adj, " (", t.test$p.signif, ")")
        }

        if (returnData)
        {
            return(a)
        }

        p1 <- ggviolin(a,
                       x        = accordingTo,
                       y        = "mean",
                       palette  = pal,
                       fill     = accordingTo,
                       color    = "mean",
                       add      = c("boxplot"),
                       outlier.shape=NA,
                       outlier.colour  = NA)
        p1 <- p1 + 
            xlab(xlab)
    }
    
    p1 <- p1 +
        stat_summary(fun=mean, #aes(color=mean),
                     geom="point", shape=20, size=4, color="red", fill="red") #+
        scale_colour_manual(label = "mean", values=c('red'),
                            name='')
#        stat_summary(aes(shape="mean"), fun.y=mean,
#                     geom="point", shape=20, size=4, color="red", fill="red") #+
#        scale_shape_manual("", values=c('Mean'='x'))
    
    pos <- position_jitter(0.2, seed=2)

    if (!is.null(y1)){
        p1 <- p1 + scale_y_continuous(limits=c(y1, y2), breaks=seq(y1, y2, by=1))
    }
    if (!is.null(compareControls)){
        require(ggrepel)
        p1 <- p1 +
            geom_hline(yintercept=meanControls + sdControls) +
            geom_hline(yintercept=meanControls - sdControls) #+
            #geom_text_repel(data = a[!a$Sample_Group == "Control" &
            #                         (a$mean > meanControls + sdControls  |
            #                          a$mean < meanControls - sdControls), ],
            #                aes(x=eval(parse(text=accordingTo)), y=mean,
            #                label=Sample),
            #                position = pos,
            #                size = 3)
        
    }
    
    if (significance){
        t.test$label <- paste0("p=", t.test$p.adj)
        print("add significance")
        p1 <- p1 + stat_pvalue_manual(t.test, label = "label",
                                      y.position = yposition[1:nrow(t.test)])
    }

    print("step12")
    p1 <- p1 +
        theme_bw() +
        ylab("M values") +
        ggtitle(title) +
#        geom_point() + 
        geom_jitter(shape=20, position=pos) + 
        theme(legend.position="right",
#              axis.text.x=element_text(size=12),
              plot.title = element_text(hjust = 0.5))
    
    return(p1)
}


data <- plot(data = data_Promoters[annotations_A2[ ! annotations_A2$Sample_Name %in% c("T0072","T0249"), ]$Sample_Name],
           Annotations = annotations_A2[ ! annotations_A2$Sample_Name %in% c("T0072","T0249"), ],
           Gene = "ATM",
           accordingTo = "Sample_Group",
           significance = TRUE,
           compareControls = TRUE,
           xlab      = "",
           y1 = -6,
           y2 = -1,
           returnData = TRUE)

meanControl = mean(data[data["Sample_Group"] == "Control",][["mean"]])
sdControl = sd(data[data["Sample_Group"] == "Control",][["mean"]])

up = meanControl + sdControl
down = meanControl - sdControl

nControls = length(data[data["Sample_Group"] == "Control",][["mean"]])
nControlsUP = sum(data[data["Sample_Group"] == "Control",][["mean"]] > up)



library(patchwork)

p1 <- plot(data = data_Promoters[annotations_A2[ ! annotations_A2$Sample_Name %in% c("T0072","T0249"), ]$Sample_Name],
           Annotations = annotations_A2[ ! annotations_A2$Sample_Name %in% c("T0072","T0249"), ],
           Gene = "ATM",
           accordingTo = "Sample_Group",
           significance = TRUE,
           compareControls = TRUE,
           xlab      = "",
           y1 = -6,
           y2 = -1)# +
#    guides(fill="none")



p1 <- p1 + #theme_Publication() + 
    theme_Publication() +
    theme(legend.position="none",
          plot.title = element_blank(),
          axis.title.x=element_blank(),
          ) +
    scale_x_discrete(labels=c("ATM heterozygous PV carriers", "Noncarriers"))





p2 <- plot(data = data_Promoters[annotations_A2[ ! annotations_A2$Sample_Name %in% c("T0072","T0249") &
                                                 annotations_A2$ATM_LOH %in% c("No LOH", "Unknown LOH"), ]$Sample_Name],
           Annotations = annotations_A2[ ! annotations_A2$Sample_Name %in% c("T0072","T0249")&
                                         annotations_A2$ATM_LOH %in% c("No LOH", "Unknown LOH"), ],
           Gene = "ATM",
           accordingTo = "Sample_Group",
           significance = TRUE,
           compareControls = TRUE,
           xlab      = "",
           y1 = -6,
           y2 = -1)# +
#    guides(fill="none")

p2 <- p2 + theme_Publication() + 
    theme(legend.position="none",
          plot.title = element_blank(),
          axis.title.x=element_blank()) +
    scale_x_discrete(labels=c("ATM heterozygous PV carriers", "Noncarriers"))

all <- p1 + p2 + plot_annotation(tag_levels = "A") +
    plot_layout(guides = 'collect', ncol=2)

ggsave(filename  = file.path(OutFilesF, "Figure 1.png"),
       plot      = all,
       dpi       = 300,
       width     =   1200,
       height    = 852,
       units     = "px",
       scale     = 3)






preComputePheatmap <- function(data,
                               annotations,
                               colAnno = c("Probable biallelic inactivation", "Sample_group"),
                               method_dist,
                               method_clust,
                               cutree_cols = 2,
                               cutree_rows = 2,
                               showRownames = F,
                               showColnames = F,
                               fontsize_row = 12)
{
    library(tidyr)
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
    
    annotations[["Probable biallelic inactivation"]] <- annotations[["Inactive"]]
    annotations$"Probable biallelic inactivation" <- factor(annotations$"Probable biallelic inactivation", levels=c("Yes", "No"))

    # Sample_Group
    annotations[['Sample_Group']] <- paste0(annotations[['Sample_Group']], '_', annotations[['label']])

    annotations <- annotations %>% mutate(Sample_Group = recode(Sample_Group,
                                                                "Non.ATM_" = "Non-ATM",
                                                                "ATM_PV"="ATM PV",
                                                                "ATM_VUS"="ATM VUS"))
    annotations["Sample group"] <- annotations["Sample_Group"]
    
    # Grade
    annotations["Grade"][annotations["Grade"] == ""] <- NA
    annotations <- annotations %>% dplyr::mutate(Grade = replace_na(Grade, "Unknown"))

    # Subtype
    #print(table(annotations$subtype))
    annotations["subtype"][annotations["subtype"] == ""] <- NA
    annotations["subtype"][annotations["subtype"] == "Missing"] <- NA
        annotations["subtype"][annotations["subtype"] == "No record"] <- NA
    annotations <- annotations %>% dplyr::mutate(subtype = replace_na(subtype, "Unknown"))
#    annotations <- annotations %>% mutate("Sample_Group" = recode(Sample_Group,
#                                                                  "ATM PV"="ATM_PV",
#                                                                  "ATM VUS"="ATM_VUS"))
    print(table(annotations$subtype))
    annoCol <- list(
        "Sample group" = c("ATM PV"="#00DAE0", "ATM VUS"="#FDD513",
                           "Non-ATM"="#E199FF"),
        "Probable biallelic inactivation" = c(Yes="#FF9289", No="#96CA00"),
        "Grade" = c("I"="yellow", "II"="orange", "III"="red", "Unknown"="white"),
        "subtype" = c("Luminal"="light blue", "Luminal A"="blue", "Luminal B"="dark blue",
                      "Luminal B/HER2+"="purple", "Triple negative"="black", "HER2+"="green",
                      "Unknown"="white")
    )
    print("italic")
    newnames <- lapply(
          rownames(data),
          function(x) bquote(italic(.(x))))    

    colAnno[which(colAnno == "Sample_Group")] <- "Sample group"
    
    heatmap <- pheatmap::pheatmap(data,
                                  treeheight_row     = 0,
                                  annotation_col     = annotations[colAnno],
                                  annotation_names_col = F,
                                  show_colnames      = showColnames,
                                  show_rownames      = showRownames,
                                  annotation_colors =  annoCol,
#                       custering_distance_cols= "correlation",
#                       custering_distance_rows= "correlation",
                                  scale              = "row",
                                  cluster_rows       = hc.features,
                                  cluster_cols       = hc.samples,
                                  cutree_cols        = cutree_cols,
                                  cutree_rows        = cutree_rows,
                                  color              = color,
                                  labels_row = as.expression(newnames),
                                  fontsize_row = fontsize_row,
                                  silent = F)

    return(heatmap)
}


###################################################################
########################## Figure 5 ###############################
###################################################################

# Vizualisation of promoters identified by the three ML algorithms
ML <- read_xlsx("/data/users/nviart/test_sel_ML_method1.xlsx")


LR <- preComputePheatmap(data=data_Promoters[na.omit(ML$'Logistic Regression'),
                                             annotations_A0A2$Sample_Name],
                         annotations   = annotations_A0A2,
                         method_dist   = "pearson",
                         method_clust  = "ward.D2",
                         cutree_cols   = 2,
                         cutree_rows   = 2,
                         showRownames  = T,
                         colAnno       = c("Probable biallelic inactivation", "Sample_Group")
                         )

print(merge(x = data.frame(cutree(LR$tree_col, k=2)), y = annotations_A0A2["Sample_Group"], by="row.names"))


RF <- preComputePheatmap(data=data_Promoters[na.omit(ML$'Random Forest'),
                                             annotations_A0A2$Sample_Name],
                         annotations   = annotations_A0A2,
                         method_dist   = "pearson",
                         method_clust  = "ward.D2",
                         cutree_cols   = 2,
                         cutree_rows   = 2,
                         showRownames  = T,
                         colAnno       = c("Probable biallelic inactivation", "Sample_Group")
                         )

XG <- preComputePheatmap(data=data_Promoters[na.omit(ML$'XGBoost'),
                                             annotations_A0A2$Sample_Name],
                         annotations   = annotations_A0A2,
                         method_dist   = "pearson",
                         method_clust  = "ward.D2",
                         cutree_cols   = 2,
                         cutree_rows   = 2,
                         showRownames  = T,
                         colAnno       = c("Probable biallelic inactivation", "Sample_Group")
                         )

# Combine plots
plot.sel.features <- patchwork::wrap_plots(list(ggplotify::as.ggplot(LR),
                                                ggplotify::as.ggplot(RF),
                                                ggplotify::as.ggplot(XG)),
                                           ncol = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = 'A')

ggsave(filename  = file.path(OutFilesF, "Figure 5.png"),
       plot      = plot.sel.features,
       dpi       = 300,
       scale     = 1,
       units     = "in",
       width     = 13,
       height    = 14)








###################################################################
########################## Figure 3 ###############################
###################################################################

source("~/ATM_Analysis/svn/analyse/script/methylation/MethPipeline/4.Interpretation.R")

# all promoters DF in A2 with all samples
H1.2 <- preComputePheatmap(data=data_Promoters[pro.A2$Gene,
                                      annotations_A0A2$Sample_Name],
                   annotations = annotations_A0A2,
                   method_dist = "pearson",
                   method_clust = "ward.D2",
                   cutree_cols        = 2,
                   cutree_rows        = 2,
                   colAnno = c("Probable biallelic inactivation", "Sample_Group")
                   )

ggsave(filename  = file.path(OutFilesF, "Figure 3.png"), 
       plot      = H1.2,
       width     = 12, #22.4
       height    = 13, #11.725,
       dpi       = 300,
       units     = "in")

