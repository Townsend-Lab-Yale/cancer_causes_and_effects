#!/bin/bash
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -J ALL_selection
#SBATCH --partition=pi_townsend
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu=15000
#SBATCH --mail-user=vincent.cannataro@yale.edu
#SBATCH --mail-type=FAIL
module load GCC/10.2.0
cd LUSC; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../selection_and_pop_level_cluster_submitter_bootstrap.R LUSC; cd ../
