#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -J fig3_data_gather
#SBATCH --partition=pi_townsend
#SBATCH -t 2:00:00
#SBATCH --mem-per-cpu=35000
#SBATCH --mail-user=vincent.cannataro@yale.edu
#SBATCH --mail-type=FAIL
module load GCC/10.2.0
~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript fig3_attributions_gather.R