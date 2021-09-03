#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

inputs <- matrix(nrow=length(args)/2,ncol=1,data=args[seq(2,length(args),2)])
rownames(inputs) <- args[seq(1,length(args),2)]
inputs[,1] <- trimws(inputs[,1])


# Create batch job and submit to the cluster
file.create(paste("selection_job_",inputs["-tumor_name",],".pbs",sep=""))
sink(file = paste("selection_job_",inputs["-tumor_name",],".pbs",sep=""))
cat("#!/bin/bash\n")
cat(paste("#SBATCH -n 1\n#SBATCH -c ",inputs["-cores",],"\n#SBATCH -J ",inputs["-tumor_name",],"_selection\n#SBATCH --partition=",inputs["-partition",],"\n#SBATCH -C haswell\n#SBATCH -t ",inputs["-time",],"\n#SBATCH --mem-per-cpu=",inputs["-mem",],"\n#SBATCH --mail-user=vincent.cannataro@yale.edu\n#SBATCH --mail-type=FAIL\n",sep=""),append=T)
cat(paste("module load libsndfile/1.0.28-GCCcore-10.2.0\n"),append = T)
cat(paste("module load ICU/67.1-GCCcore-10.2.0\n"),append = T)
cat(paste("~/../../pi/townsend/vlc24/R_source/R-4.0.2/bin/Rscript ../selection_and_pop_level_cluster_mp.R",paste(args,collapse = " "),sep=" "),append = T)
sink()
system(paste("sbatch selection_job_",inputs["-tumor_name",],".pbs",sep=""))
