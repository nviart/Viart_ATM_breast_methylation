
# Set the path to renv project
renv.path = "/data/users/nviart/ATM_Analysis/svn/analyse/script/"

# Set the path to where you want to save the results of the scripts
Origin = "/data/users/nviart/ATM_Analysis/svn/analyse/results/FinalMethylation"

# Set the path to the Illumina manifest
manifest.path = "/data/users/nviart/ATM_Analysis/data/methylation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"

# Set to TRUE if you only want to use samples for which methylation was measured with the EPIC array
onlyEPIC = FALSE

# Set the paths to where are located the methylation data and the descriptive csv files 
EpicFolder = "/data/users/nviart/ATM_Analysis/data/methylation"
DescriptorFileName1 = "Annotations_30.11.23_Curie.csv"
DescriptorFileName2 = "Annotations_13.02.2024_Australia_ABCFR_MCCS.csv"


# Set to FALSE if you do not want to suppress probes from https://github.com/sirselim/illumina450k_filtering
probeFiltering = TRUE
# Set the path to where the files from https://github.com/sirselim/illumina450k_filtering are located
path.filters = "/data/users/nviart/ATM_Analysis/data/methylation/Probes_filtering/"
