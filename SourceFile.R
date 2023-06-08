#############################
# this is a template source file
# please change all paths accordingly
#############################

#############################
# Working directory
#############################
basicpath = "/path_to_your_project/"
basicpath_scripts = "/path_to_your_project/scripts/"

#############################
# R library and R packages
#############################
.libPaths("/path_to_your_R_libraries/") 
.libPaths()

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(MendelianRandomization))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(WriteXLS))

#############################
# Downloaded data sets 
#############################
path_downloaded = paste0(basicpath, "reference_data/")
