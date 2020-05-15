#############################################################################################################################
##                                                                                                                      
##  Create genomedb.list for the mtdna_processing.sh script                                                                                                         
##                                                                                                                      
##  Date: 24 July 2019                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##  Summary: genomedb.list is a tab delimited file with <sample ID> <full path to corresponding .g.vcf.gz>                                                                                                                 
##                                                                                                                      
############################################################################################################################

# clear workspace beforehand
rm(list = ls())

# load required libraries
library(tidyverse)
library(reshape2)

# get arguments following the script
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nError!\nusage: Rscript make_genomedbList.R /path/to/output_dir\n")
  quit()
}

# set output directory
odir <- args[1]     # odir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS/output"
message(odir)

# list all gvcfs
vcf_files <- list.files(paste(odir, "/gatk/gvcfs", sep = ""), pattern = ".g.vcf.gz$", full.names = T,recursive = T)

p <- str_split_fixed(basename(vcf_files), "\\.", 2)[,1]

# create data frame with column ID and BAM 
this <- data.frame(ID=p,
                   VCF=vcf_files, 
                   stringsAsFactors = F)
# write talbe
write.table(this, paste(odir, "/gatk/gvcfs/genomedb.list", sep = ""), 
            sep = "\t", col.names = F, row.names = F,quote = F)
