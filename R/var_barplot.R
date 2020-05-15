#############################################################################################################################
##                                                                                                                      
##  plot bar plots for the total number of mutations detected in a patient x - #muts; y - patients                                                                                                
##                                                                                                                      
##  Date: 04.02.2020                                                                                                              
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                                                                                                                                            
############################################################################################################################

# clear workspace beforehand
rm(list = ls())

# load required libraries
library(tidyverse)
library(reshape2)
library(stringr)
library(MutationalPatterns)

# get arguments following the script
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=4){
  message("\nError!\nusage: Rscript var_barplot.R /path/to/all_patient_variant_file\n /path/to/output_dir\n")
  quit()
}

# set output directory
odir <- args[2]     # odir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K/plots"
message(odir)
dir.create(odir)
setwd(odir)

# read in data_frame 
allvar.df <- read.table(args[1], header = T) # allvar.df <- read.table("/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K/mtdnaserver/filtered/all_patient_0.05_0.01_filtered.var.txt", header = T)

# get the thresholds
ts <- paste(str_split_fixed(basename(args[1]), "_", 5)[3], str_split_fixed(basename(args[1]), "_", 5)[4], sep = "_") # args <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K/mtdnaserver/filtered/all_patient_0.05_0.01_filtered.var.txt"

# get patient ids 
patients <- unique(str_split_fixed(colnames(allvar.df)[seq(2, ncol(allvar.df), 2)], "_", 2)[,1])
patient_list <- str_split_fixed(patients, "\\.", 2)[,2]

# count the number of mutations per patient
# patient <- patient_list[1]
n_mut.list <- c()
for (patient in patient_list){
  
  # which patient?
  message(patient)
  
  # grep all cols in df which belong to this patient
  p_cols <- grep(pattern = patient, colnames(allvar.df)) 
  
  # i <- 2
  n_muts <- 0
  for (i in seq(2, ncol(allvar.df),2)){
    if (i %in% p_cols[(length(p_cols)-1)]){
      n_muts <- n_muts ++ sum(!is.na(allvar.df[,i]))
    } else {
      next
    }
  }
  this <- data.frame(Patient=patient, 
                     Number=n_muts, 
                     stringsAsFactors = F)
  n_mut.list <- rbind(n_mut.list, this)
}

# plot patient barplot  
ggplot(n_mut.list, aes(Patient, Number))+ 
  geom_bar(aes(fill=Patient), stat = "identity") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Mutations per Patient") +
  xlab("Patients") +
  ylab("# Mutations") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 9, angle = 90, hjust = .5, vjust = .5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 9, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 12, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 15, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 15, hjust = .5, face = "bold")) +
  theme(legend.position="none")

ggsave(file.path(odir, paste("/Barplot_Mutations_Patients_", ts, ".pdf", sep = "")), width = 10, height = 10)

