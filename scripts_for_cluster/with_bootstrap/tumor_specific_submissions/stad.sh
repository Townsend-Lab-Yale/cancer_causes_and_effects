#!/bin/bash
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -J STAD_bootstrap
#SBATCH --partition=general
#SBATCH -t 10:00:00
#SBATCH --mem-per-cpu=16000
#SBATCH --mail-user=vincent.cannataro@yale.edu
#SBATCH --mail-type=FAIL
module load GCC/10.2.0
cd STAD; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../selection_and_pop_level_cluster_submitter_bootstrap.R STAD; cd ../

