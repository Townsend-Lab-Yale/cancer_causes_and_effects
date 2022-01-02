#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



job_ind <- 1

library(tidyverse)
library(cancereffectsizeR)


source("../parameters_for_tumor_types.R")

parameters_for_run <- parameters_for_run[[args]]

n_cores <- parameters_for_run %>% 
  filter(input_parameter == "n_cores") %>% 
  pull(input_value)

n_boot <- parameters_for_run %>% 
  filter(input_parameter == "n_boot") %>% 
  pull(input_value)

n_split <- parameters_for_run %>% 
  filter(input_parameter == "split") %>% 
  pull(input_value)


tumor_name <- parameters_for_run %>% 
  filter(input_parameter == "tumor_name") %>% 
  pull(input_value)


partition <- parameters_for_run %>% 
  filter(input_parameter == "partition") %>% 
  pull(input_value)


time <- parameters_for_run %>% 
  filter(input_parameter == "time") %>% 
  pull(input_value)

memory <- parameters_for_run %>% 
  filter(input_parameter == "memory") %>% 
  pull(input_value)


file.create(paste("combine_job_",args,"_",job_ind,".sh",sep=""))
sink(file = paste("combine_job_",args,"_",job_ind,".sh",sep=""))
cat("#!/bin/bash\n")
cat(paste("#SBATCH -n 1\n#SBATCH -c ",n_cores,"\n#SBATCH -J ", tumor_name ,"_selection\n#SBATCH --partition=",partition,"\n#SBATCH -t ", time ,"\n#SBATCH --mem-per-cpu=",memory,"\n#SBATCH --mail-user=vincent.cannataro@yale.edu\n#SBATCH --mail-type=FAIL\n",sep=""),append=T)
cat(paste("module load GCC/10.2.0\n"),append = T)
cat(paste("~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner.R",args,job_ind,sep=" "),append = T)
sink()
system(paste("sbatch combine_job_",args,"_",job_ind,".sh",sep=""))
