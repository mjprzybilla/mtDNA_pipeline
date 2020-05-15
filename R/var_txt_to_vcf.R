#############################################################################################################################
##                                                                                                                      
##  Convert the filtered mtdnaserver output into a vcf file                                                                                                          
##                                                                                                                      
##  Date: 16 August 2019                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##  Summary: This script generates a vcf file based on the filtered mutation file from mtdnaserver.
##           The file is used to run mode 1 of cellSNP for genotyping sc data.  
##                                                                                                                              
##                                                                                                                      
############################################################################################################################

# clear workspace beforehand
rm(list = ls())

# load required libraries
library(tidyverse)
library(reshape2)
library(ggplot2)

# get arguments following the script
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=4){
  message("\nError!\nusage: Rscript var_txt_to_vcf.R /path/to/patient_variant_files\n /path/to/output_dir\n control_threshold\n tumor_threshold\n")
  quit()
}

# Rscript /home/przybilm/bsub_command/WGS/mtdna_pipeline/var_txt_to_vcf.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_067/mtdnaserver/filtered /icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_067/mtdnaserver/filtered/vcf 0.05 0.01

# get input directory containing the output from mtdnaserver
wdir <- args[1]     # wdir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K/mtdnaserver/filtered"
# wdir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS/output/mtdnaserver/filtered"

# set output directory
odir <- args[2]     # odir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K/mtdnaserver/filtered/vcf"
# odir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS/output/mtdnaserver/filtered/vcf"
message(odir)
dir.create(odir)
setwd(odir)

# set thresholds for controls and tumours 
ctrl_t <- args[3]   # ctrl_t <- 0.5
tumor_t <- args[4]  # tumor_t <- 0.05

# list all files containing variants called with mtdnaserver
allvarfiles <- list.files(wdir, pattern = paste("_", ctrl_t, "_", tumor_t, "_filtered.var.txt", sep = ""), full.names = T)
# allvarfiles <- allvarfiles[-grep("all_patient", allvarfiles)]

# create vcf file for cellSNP input 
# file <- allvarfiles[1]
for (file in allvarfiles){
  
  # message the file we are dealing with
  message(file)
  
  # read it in
  table.tmp <- read.table(file, header = T)
  
  # get patient id 
  patient <- str_split_fixed(basename(file), "_", 2)[1]
  
  # read in REF mutations from https://www.mitomap.org/foswiki/bin/view/MITOMAP/MitoSeqs
  ref.vcf <- read.delim("~/mtdna_pipeline/REFs/chrM_REF_mutations.vcf", header = T)
  colnames(ref.vcf)[1] <- "#CHROM"
  
  # iterate over each row and convert it to a vcf file
  vcf.list <- c()
  if (nrow(table.tmp) > 0){
    for (i in 1:nrow(table.tmp)){
      split <- str_split_fixed(table.tmp[i,1], ":", 3)
      vcf <- data.frame(CHROM = "chrM", 
                        POS = as.numeric(split[1]), 
                        ID = ".",
                        REF = split[2],
                        ALT = split[3],
                        QUAL = ".",
                        FILTER = ".", 
                        INFO = ".", stringsAsFactors = F)
      colnames(vcf)[1] <- "#CHROM"
      vcf.list <- rbind(vcf.list, vcf)
    }
    
    # negate the %in% operator
    `%notin%` <- Negate(`%in%`)
    
    # add a filtering steps of variants
    vcf.list <- vcf.list[vcf.list$POS %notin% ref.vcf$POS,]
    
    # write vcf list for each individual patient to file
    write.table(vcf.list, paste(odir, "/", patient, "_", ctrl_t, "_", tumor_t, "_filtered.var.vcf", sep = ""), row.names = F, quote = F, sep = "\t")
  }
}

