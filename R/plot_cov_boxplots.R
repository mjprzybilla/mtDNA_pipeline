#############################################################################################################################
##                                                                                                                      
##  Plotting of the coverage across the mt genome after filtering the raw mitochondrial bam files                                                                                                    
##                                                                                                                      
##  Date: 27 August 2019                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##  Summary: This script makes use of the Rtsne package to plot the coverage of the WGS bams for the mt genome.                                                           
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list = ls())

# load packages
library(ShortRead)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(tidyverse)
library(reshape2)

# get arguments following the script
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=4){
  message("\nError!\nusage: Rscript plot_mtdna_coverage.R /path/to/filtered_mt_bams\n /path/to/output_dir\n control_threshold\n tumor_threshold\n")
  quit()
}

# set directory
wdir <- args[1]   # wdir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS/output/filtered_unique_rg"

# set output directory
odir <- args[2]     # odir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS/output/plots"
message(odir)
dir.create(odir)
setwd(odir)

# =============================================================================================================================
# coverage boxplot (after filtering and remapping)
# =============================================================================================================================
rm(list = ls())
library(Rsamtools)
library(GenomicAlignments)

# construct a list of mtDNA bam files
filtered_file_names <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/przybilm/filtered_unique_rg", full.names = T, recursive = T, pattern = ".bam$")
filtered_K08K_names <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K/filtered_unique_rg", full.names = T, recursive = T, pattern = ".bam$")
total_filtered_names <- c(filtered_file_names, filtered_K08K_names)

mtdna_df <- read.table("/icgc/dkfzlsdf/analysis/B260/projects/przybilm/combined_mtdnaserver/mtdnaserver_all.var.txt", header = T)

# patients 
patients <- strsplit(colnames(mtdna_df)[seq(3, length(colnames(mtdna_df)),2)], "_")
patient_list <- c()
#i <- 1
for (i in 1:length(patients)){
  patient <- patients[[i]][1]
  patient <- gsub("\\.", "-", patient)
  patient_list <- rbind(patient_list, patient)
}
patient_list <- unique(patient_list)

ATAC_file_names <- c()
# patient <- patient_list[1]
for (patient in patient_list){
  files <- grep(patient, total_filtered_names, value = T)
  ATAC_file_names <- c(ATAC_file_names, files)
}

coverage_list <- list()
# loop over list of mtDNA bam files
for (bam in filtered_file_names){
  message(bam)
  chrid=basename(bam)
  
  if (grepl("*control*", chrid)) {
    type <- "control"
  } else {
    type <- "tumor"
  }
  
  # Read alignments
  bm=readGAlignments(bam)
  
  # Calculate coverage at each position on each chromosome
  dat=coverage(bm)
  
  # Take the reads for chromosome 11 and thin (otherwise plots will take a lot of memory)
  chrname="MT"
  rawreads=as.numeric(dat[[chrname]])
  # determine the segments, intervals for the plotting
  inds=seq(1,length(rawreads))
  averagereads = mean(rawreads)
  
  # initialize dataframe with ID and type
  this <- data.frame(ID=chrid,
                     TYPE=type,
                     stringsAsFactors = F)
  
  # dataframe with reads per position in bp
  df <- as.data.frame(rbind(inds, rawreads))
  colnames(df) = df[1, ] # the first row will be the header
  df = df[-1, ]      
  
  # merge both dataframes to have one row per sample 
  patient_row <- merge(this, df)
  patient_row <- merge(patient_row, averagereads)
  colnames(patient_row)[16572] <- "MEAN"
  
  # bind this together with coverage_list
  coverage_list <- rbind(coverage_list,patient_row)
}

cov_df <- coverage_list[, c(1,2,16572)]

tumor <- count(coverage_list$TYPE == "tumor")
controls <- count(coverage_list$TYPE == "control")

# box plot 
g <- ggplot(coverage_list, aes(TYPE, MEAN)) +
  geom_boxplot(aes(fill=factor(TYPE))) + 
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .4, 
               fill="black") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Coverage box plot - Filtered bams",
       x="Type",
       y="Average coverage") + 
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 0, hjust = .5, vjust = .5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 7, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 12, hjust = .5, face = "bold"), 
        legend.position = "none") +
  scale_fill_brewer(palette="Dark2")

ggsave(file.path("/home/przybilm/mtDNA_plots/filtered_unique_rg/", paste("Filtered_Coverage_boxplot", ".pdf", sep = "")))

# =============================================================================================================================
# coverage boxplot of the ATAC-seq patients
# =============================================================================================================================

ATAC_coverage_list <- list()
# loop over list of mtDNA bam files
for (bam in ATAC_file_names){
  message(bam)
  chrid=basename(bam)
  
  if (grepl("*control*", chrid)) {
    type <- "control"
  } else {
    type <- "tumor"
  }
  
  # Read alignments
  bm=readGAlignments(bam)
  
  # Calculate coverage at each position on each chromosome
  dat=coverage(bm)
  
  # Take the reads for chromosome 11 and thin (otherwise plots will take a lot of memory)
  chrname="MT"
  rawreads=as.numeric(dat[[chrname]])
  # determine the segments, intervals for the plotting
  inds=seq(1,length(rawreads))
  averagereads = mean(rawreads)
  
  # initialize dataframe with ID and type
  this <- data.frame(ID=chrid,
                     TYPE=type,
                     stringsAsFactors = F)
  
  # dataframe with reads per position in bp
  df <- as.data.frame(rbind(inds, rawreads))
  colnames(df) = df[1, ] # the first row will be the header
  df = df[-1, ]      
  
  # merge both dataframes to have one row per sample 
  patient_row <- merge(this, df)
  patient_row <- merge(patient_row, averagereads)
  colnames(patient_row)[16572] <- "MEAN"
  
  # bind this together with coverage_list
  ATAC_coverage_list <- rbind(ATAC_coverage_list,patient_row)
}

grep("control", ATAC_coverage_list$ID, value = T)

tumor <- count(ATAC_coverage_list$TYPE == "tumor")
controls <- count(ATAC_coverage_list$TYPE == "control")

# box plot 
g_ATAC <- ggplot(ATAC_coverage_list, aes(TYPE, MEAN)) +
  geom_boxplot(aes(fill=factor(TYPE))) + 
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .4, 
               fill="black") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="Type",
       y="Mean coverage") + 
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 0, hjust = .5, vjust = .5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 7, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 12, hjust = .5, face = "bold"), 
        legend.position = "none") +
  scale_fill_brewer(palette="Dark2")

ggsave("/home/przybilm/mtDNA_plots/filtered_unique_rg/ATAC_patients_coverageBoxplot.pdf", height = 10, width = 10)
