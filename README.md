#  Pre-processing pipeline of methylation data (Illumina 450K and EPIC arrays)

This repository contains the code used for the paper Viart *et al.* (2025) entitled "Breast tumors from ATM pathogenic variant carriers display a specific genome-wide DNA methylation profile"

## 1- Clone the repository
To be able to use the code of this repository, you have first to clone the repository with the command:
```
git clone https://github.com/nviart/Viart_ATM_breast_methylation.git
```

## 2- Package installation
This code is mainly written in R language. The easiest way to use it with the same development environment is to use the versionning system renv (https://rstudio.github.io/renv/index.html). First install renv package if not done. Then, change in the config.R file the path to renv environment ("renv.path" variable). It will allow to load the environment in each of the other scripts.

## Download the recquired files to run the pipeline


## 2- Modify the config file
The config.R file contains some variables that will be used by the other scripts. Set them according to your own environment.
You have to indicate:
* `onlyEPIC`: if only EPIC arrays will be used
* `EPICfolder`: the folder containig ALL your methylation arrays data
* `DescriptorFileName1` and `DescriptorFileName1`: the names of the csv descriptive files of the samples that should contain at least ....
* `probeFiltering`: boolean, if you want to filter the probes described in Pidsley *et al.* (2016) 
* `path.filters`: path to the files downloaded from https://github.com/sirselim/illumina450k_filtering


function.path = "/data/users/nviart/ATM_Analysis/svn/analyse/script/methylation/MethPipeline/"~