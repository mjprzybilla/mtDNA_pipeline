#############################################################################################################################
##                                                                                                                      
##  Filter the combined mtDNA server output                                                                                                         
##                                                                                                                      
##  Date: 16 August 2019                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##  Summary: This script applys filter to the bulk WGS data from mtdnaserver
##                                                                                                                              
##                                                                                                                      
############################################################################################################################

# clear workspace
rm(list = ls())

# load libraries
library(stringr)
library(ggplot2)
library(tidyverse)
library(reshape2)

# get arguments following the script
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=4){
  message("\nError!\nusage: Rscript mtdna_var_filtering.R /path/to/mtdnaserver_output\n /path/to/output_dir\n control_threshold\n tumor_threshold\n")
  quit()
}

# get input directory containing the output from mtdnaserver
wdir <- args[1]     # wdir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS/output/mtdnaserver"

# set output directory
odir <- args[2]     # odir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS/output/mtdnaserver/filtered"
message(odir)
dir.create(odir)
setwd(odir)

# list all files containing variants called with mtdnaserver
allvarfiles <- list.files(wdir, pattern = "_all.var.txt$", full.names = T)
allvarfiles <- allvarfiles[-grep("all_patient", allvarfiles)]

# set a first threshold for both VAFs and a second one in case there is only one sample
t1 <- args[3]   # t1 <- 0.5
t2 <- args[4]  # t2 <- 0.05

# file <- allvarfiles[1]
for (file in allvarfiles){
  
  # message file of interest
  message(file)
  
  # read in patient file and replace NAs with 0s
  table.tmp <- read.table(file, header = T)
  table.tmp[is.na(table.tmp)] <- 0
  
  # grep columns that belong to a tumor or a control sample
  sample1_col <- grep(pattern = "Rep1", colnames(table.tmp),value = T)
  sample2_col <- grep(pattern = "Rep2", colnames(table.tmp),value = T)
  
  # get information about the VAF and the corresponding coverage at the position
  vcf_sample1_col <- grep(sample1_col,pattern = "_VariantLevel",value = T)
  vcf_sample2_col <- grep(sample2_col,pattern = "_VariantLevel",value = T)
  cov_sample1_col <- grep(sample1_col,pattern = "_Coverage",value = T)
  cov_sample2_col <- grep(sample2_col,pattern = "_Coverage",value = T)
  
  # get patient name
  patient <- paste(str_split_fixed(basename(file), "_", 3)[,1], str_split_fixed(basename(file), "_", 3)[,2], sep = "_")
  
  # i = 1
  keep <- c()
  for(i in 1:nrow(table.tmp)){
    
    # look at one point mutation at a time
    pm <- table.tmp[i,]
    
    # get the variant allele frequency of the mutation to apply the filter
    vafs_sample1 <- as.numeric(pm[,vcf_sample1_col])
    vafs_sample2 <- as.numeric(pm[,vcf_sample2_col])
    
    # check whether the variants exceed the thresholds
    if(all(vafs_sample1 <= t1 & vafs_sample2 <= t1)){
      if(length(which(vafs_sample1 >= t2 | vafs_sample2 >= t2)) >= 1){
        keep <- c(keep, i) 
      }
    }
  }
  
  # subset the original table by the mutations to keep
  table_filter <- table.tmp[keep, ]
  
  # write to file
  write.table(table_filter, paste(odir, "/" , patient, "_", t1, "_", t2, "_filtered.var.txt", sep = ""), 
              sep = "\t", row.names = F, quote = F)
}

#####################################################################################
# merge all filtered patient-specific variants together
#####################################################################################

# list all files filtered with the identical thresholds
filtered_allvarfiles <- list.files(odir, pattern = paste(ctrl_t, "_", tumor_t, "_filtered.var.txt$", sep = ""), full.names = T)

# read files in
filtered_allvar.dfs <- lapply(filtered_allvarfiles, read.table, header = T)

# merge together and write to an individual file
filtered_allvar.dfs <- Reduce(function(x, y) merge(x, y, by = "Mutation", all=TRUE), filtered_allvar.dfs)
write.table(filtered_allvar.dfs, paste(odir, "/all_patient_", ctrl_t, "_", tumor_t, "_filtered.var.txt", sep = ""), row.names = F, quote = F, sep = "\t")
