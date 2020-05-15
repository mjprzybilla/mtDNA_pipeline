#############################################################################################################################
##                                                                                                                      
##  Create a samplesheet with the name of the sample and its respective file path as input for the 
##  mtdna_preprocessing_variantcalling.sh script                                                                                                      
##                                                                                                                      
##  Date: 30.01.2020                                                                                                           
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                                                                                                                                        
############################################################################################################################

# clear workspace beforehand
rm(list = ls())

# load required libraries
library(tidyverse)
library(reshape2)

# get arguments following the script
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2){
  message("\nError!\nusage: Rscript get_sampleSheet.R /path/to/bamfile_directory\n /path/to/output_dir\n")
  quit()
}

# get input directory
wdir <- args[1]     # wdir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS"

# set output directory
odir <- args[2]     # odir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS/output"
message(odir)
setwd(odir)

# list all patient folders with ids
patients <- unique(list.files(wdir,full.names = T))

# create vector bamfile_list
bamfile_list <- c()

# loop over ids (patients) in list files 
# p <- patients[1]
for(p in patients){
  message(p)
  # grep the directories which are control or tumor
  samples <- grep(list.files(p,full.names = T),pattern = "control|tumor|Rep",value = T)
  # loop over the control or tumor directories
  # x <- samples[1]
  for(x in samples){
    message(x)
    # list the bam files for the respective control or tumour sample
    bam <- list.files(x,pattern = "\\.bam$",full.names = T,recursive = T)
    # create data frame with column ID and BAM 
    this <- data.frame(ID=paste(basename(p),basename(x),sep = "_"),
                       BAM=bam,
                       stringsAsFactors = F)
    # bind this together with bamfile_list
    bamfile_list <- rbind(bamfile_list,this)
  }
}

# check for duplicates
# dups <- bamfile_list$ID[which(duplicated(bamfile_list$ID))]
# bamfile_list <- bamfile_list[which(bamfile_list$ID != dups[1]),]

# write bamfile_list as txt file
write.table(bamfile_list, "bamlist.txt", sep = "\t", col.names = F, row.names = F,quote = F)
