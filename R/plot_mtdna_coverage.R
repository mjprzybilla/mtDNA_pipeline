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

# construct a list of mtDNA bam files
filtered_mt.bam <- list.files(wdir, full.names = T, recursive = T, pattern = ".bam$")

# get patient ids
patients <- unique(paste(str_split_fixed(basename(filtered_mt.bam), "_", 3)[,1], str_split_fixed(basename(filtered_mt.bam), "_", 3)[,2],  sep = "_"))

# create vector with individual bases of the mt genome
cov.list <- c()
inds=seq(1,16569)
cov.list <- cbind(cov.list, inds)
colnames(cov.list)[1] <- "Position"

# loop over list of mtDNA bam files
# bam.tmp <- filtered_mt.bam[1]
for (bam.tmp in filtered_file_names){
  
  # which bam file?
  message(bam.tmp)
  
  # get chrid for file
  file.id <- str_split_fixed(basename(bam.tmp), "\\.", 2)[1]
  
  # check if its control tumor or normal Rep
  if (grepl("*control*", file.id)) {
    type <- "control"
  } else if (grepl("tumor", file.id)) {
    type <- "tumor"
  } else {
    type <- "normal"
  }
  
  # Read alignments
  bm=readGAlignments(bam.tmp)
  
  # Calculate coverage at each position on each chromosome
  dat=coverage(bm)
  
  # Take the reads for the mitochondrial genome
  chrname="MT"
  rawreads=as.numeric(dat[[chrname]])
  
  # determine the segments, intervals for the plotting
  average.reads = mean(rawreads)
  
  # initialize dataframe with ID and type
  this <- data.frame(file.id = rawreads,
                     stringsAsFactors = F)
  colnames(this) <- file.id
  
  # dataframe with reads per position in bp
  cov.list <- cbind(cov.list, this)
}

# calculate the median coverage over all positions
cov.median <- colMedians(as.matrix(cov.list[,c(2:ncol(cov.list))]))
names(cov.median) <- colnames(cov.list)[c(2:ncol(cov.list))]

# calculate the average coverage for each patient
cov.mean <- colMeans(as.matrix(cov.list[,c(2:ncol(cov.list))]))
names(cov.mean) <- colnames(cov.list)[c(2:ncol(cov.list))]

# plot the coverage for each patient 
# patient <- patients[1]
for (patient in patients){
  
  # which patient?
  message(patient)
  
  # number of columns per patient (aka samples) with the respective sample names
  patient_col <- grep(pattern = patient, colnames(cov.list))
  name_col <- grep(pattern = patient, colnames(cov.list), value = T)
  n_pc <- length(patient_col)
  
  # plot depending on the number of samples
  if(n_pc == 2){
    Y1 <- cov.list[, patient_col[1]]
    Y2 <- cov.list[, patient_col[2]]
      
    g <- ggplot(cov.list, aes(x=cov.list$Position)) +
      geom_line(aes(y= Y1)) +
      geom_line(aes(y= Y2), colour ="red") +
      labs(title = patient, 
           y = "Coverage",
           x = "mtDNA [bp]", 
           fill = "Samples") +
      theme_classic() +
      theme(axis.text.x = element_text(vjust = 0.5)) +
      # scale_y_continuous(limits= c(0, 100000), breaks = seq(0, 100000, by = 5000)) + 
      scale_x_continuous(breaks = seq(0, 16569, by = 1000)) +
      guides(fill = name_col) +
      # geom_hline(yintercept=median, color = "grey", linetype = "dashed") +
      # geom_text(aes(0,median), label = paste("Median coverage =", round(median,2), sep= " "), vjust = -10, hjust = -1) + 
      theme(axis.text.x = element_text(color = "black", size = 9, angle = 0, hjust = .5, vjust = .5, face = "bold"),
            axis.text.y = element_text(color = "black", size = 9, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
            axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "bold"),
            axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                        margin = margin(t = 0, r = 15, b = 0, l = 0)), 
            plot.title = element_text(color = "black", size = 20, hjust = .5, face = "bold"))
    
    # print plot
    print(g)
    
  } else {
    Y1 <- cov.list[,patient_col[1]]
    Y2 <- cov.list[,patient_col[2]]
    Y3 <- cov.list[,patient_col[3]]
      
    g <- ggplot(cov.list, aes(x=cov.list$Position)) +
      geom_line(aes(y= Y1))+
      geom_line(aes(y= Y2), colour ="red") +
      geom_line(aes(y= Y3), colour ="blue") +
      labs(title = patient[[1]], 
           y = "Coverage",
           x = "mtDNA [bp]",
           fill = "Samples") +
      theme_classic() +
      theme(axis.text.x = element_text(vjust = 0.5)) +
      # scale_y_continuous(limits= c(0, 100000), breaks = seq(0, 100000, by = 5000)) + 
      scale_x_continuous(breaks = seq(0, 16569, by = 1000)) +
      # geom_hline(yintercept=median, color = "grey", linetype = "dashed") +
      # geom_text(aes(0,median), label = paste("Median coverage =", round(median,2), sep= " "), vjust = -10, hjust = -1) + 
      theme(axis.text.x = element_text(color = "black", size = 9, angle = 0, hjust = .5, vjust = .5, face = "bold"),
            axis.text.y = element_text(color = "black", size = 9, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
            axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "bold"),
            axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                        margin = margin(t = 0, r = 15, b = 0, l = 0)), 
            plot.title = element_text(color = "black", size = 20, hjust = .5, face = "bold")) +
      guide_legend(name_col)
    
    # print plot
    print(g)
    
  }
  
  ggsave(paste(odir, "/", patient, "_multi_coverage.pdf", sep = ""), width = 11, height = 11)
  
}
