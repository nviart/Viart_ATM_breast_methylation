
# Set the path to renv project
renv.path = "/data/users/nviart/ATM_Analysis/svn/analyse/script/"

# Set the path to where you want to save the results of the scripts
Origin = "/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation"
#FinalMethylation_onlyEPIC"
# FinalMethylation_withoutProbesFiltering"

# Set to TRUE if you only want to use samples for which methylation was measured with the EPIC array
onlyEPIC = FALSE

# Set the paths to where are located the methylation data and the descriptive csv files 
EpicFolder = "/data/users/nviart/ATM_Analysis/data/methylation"
DescriptorFileName1 = "Annotations_30.11.23_Curie.csv"
DescriptorFileName2 = "Annotations_13.02.2024_Australia_ABCFR_MCCS.csv"


# Set to FALSE if you do not want to suppress probes from https://github.com/sirselim/illumina450k_filtering
probeFiltering = TRUE
# Set the path to where the files from https://github.com/sirselim/illumina450k_filtering are located
path.filters = 

function.path = "/data/users/nviart/ATM_Analysis/svn/analyse/script/methylation/MethPipeline/"




# For TCGA config
TCGA = TRUE
Origin = "/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation_TCGA_balanced/"
# Set to TRUE if you only want to use samples for which methylation was measured with the EPIC array
onlyEPIC = FALSE

# Set to FALSE if you do not want to suppress probes from https://github.com/sirselim/illumina450k_filtering
probeFiltering = TRUE

renv.path = "/data/users/nviart/ATM_Analysis/svn/analyse/script/"






#scriptPath <- function() {
#    cmdArgs <- commandArgs(trailingOnly = FALSE)
#    needle <- "--file="
#    match <- grep(needle, cmdArgs)
#    if (length(match) > 0) {
#                                        # Rscript
#        return(normalizePath(sub(needle, "", cmdArgs[match])))
#    } else {
#                                        # 'source'd via R console
#        return(normalizePath(sys.frames()[[1]]$ofile))
#    }
#}

#scriptFolder <- dirname(scriptPath())
