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

# Rscript /home/przybilm/bsub_command/WGS/mtdna_pipeline/mtdna_var_filtering.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_067/mtdnaserver /icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_067/mtdnaserver/filtered 0.05 0.01

# get input directory containing the output from mtdnaserver
wdir <- args[1]     # wdir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K/mtdnaserver/lowVAF"

# set output directory
odir <- args[2]     # odir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K/mtdnaserver/filtered"
message(odir)
dir.create(odir)
setwd(odir)

# list all files containing variants called with mtdnaserver
allvarfiles <- list.files(wdir, pattern = "_all.var.txt$", full.names = T)
# allvarfiles <- allvarfiles[-grep("all_patient", allvarfiles)]

# set thresholds for controls and tumours 
ctrl_t <- args[3]   # ctrl_t <- 0.05
message(ctrl_t)
tumor_t <- args[4]  # tumor_t <- 0.01
message(tumor_t)

# file <- allvarfiles[1]
for (file in allvarfiles){
  
  # message file of interest
  message(file)
  
  # read in patient file and replace NAs with 0s
  table.tmp <- read.table(file, header = T)
  table.tmp[is.na(table.tmp)] <- 0
  
  # grep columns that belong to a tumor or a control sample
  tumor_col <- grep(pattern = "tumor", colnames(table.tmp),value = T)
  control_col <- grep(pattern = "control", colnames(table.tmp),value = T)
  
  # get information about the VAF and the corresponding coverage at the position
  vcf_ctrl_col <- grep(control_col,pattern = "_VariantLevel",value = T)
  vcf_tumor_col <- grep(tumor_col,pattern = "_VariantLevel",value = T)
  cov_ctrl_col <- grep(control_col,pattern = "_Coverage",value = T)
  cov_tumor_col <- grep(tumor_col,pattern = "_Coverage",value = T)
  
  # get patient name
  patient <- str_split_fixed(basename(file), "_", 2)[,1]
  
  # i = 1
  keep <- c()
  for(i in 1:nrow(table.tmp)){
    
    # look at one point mutation at a time
    pm <- table.tmp[i,]
    
    # get the variant allele frequency of the mutation to apply the filter
    vafs_ctrl <- as.numeric(pm[,vcf_ctrl_col])
    vafs_tumor <- as.numeric(pm[,vcf_tumor_col])
    
    # check whether the variants exceed the thresholds
    if(all(vafs_ctrl <= ctrl_t)){
      if(length(which(vafs_tumor >= tumor_t)) >= 1){
        keep <- c(keep, i) 
      }
    }
  }
  
  # subset the original table by the mutations to keep
  table_filter <- table.tmp[keep, ]
  
  # write to file
  write.table(table_filter, paste(odir, "/" , patient,"_", ctrl_t, "_", tumor_t, "_filtered.var.txt", sep = ""), 
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
