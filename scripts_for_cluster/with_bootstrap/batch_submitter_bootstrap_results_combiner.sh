#!/bin/bash
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -J ALL_selection_combine
#SBATCH --partition=pi_townsend
#SBATCH -t 10:00:00
#SBATCH --mem-per-cpu=16000
#SBATCH --mail-user=vincent.cannataro@yale.edu
#SBATCH --mail-type=FAIL
module load GCC/10.2.0
cd UCEC; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R UCEC_bootcombine; cd ../
cd LUSC; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R LUSC; cd ../
cd BLCA; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R BLCA; cd ../
cd BRCA_ER_neg; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R BRCA_ER_neg; cd ../
cd BRCA_ER_pos; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R BRCA_ER_pos; cd ../
cd CESC; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R CESC; cd ../
cd COAD; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R COAD; cd ../
cd ESCA; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R ESCA; cd ../
cd ESCC; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R ESCC; cd ../
cd GBM; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R GBM; cd ../
cd HNSC_HPVneg; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R HNSC_HPVneg; cd ../
cd HNSC_HPVpos; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R HNSC_HPVpos; cd ../
cd KIRC; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R KIRC; cd ../
cd LGG; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R LGG; cd ../
cd LUAD; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R LUAD; cd ../
cd LIHC; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R LIHC; cd ../
cd OV; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R OV; cd ../
cd PAAD; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R PAAD; cd ../
cd PRAD; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R PRAD; cd ../
cd READ; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R READ; cd ../
cd SKCM_primary; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R SKCM_primary; cd ../
cd SKCM_metastasis; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R SKCM_metastasis; cd ../
cd STAD; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R STAD; cd ../
cd THCA; ~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../bootstrap_results_combiner_submit.R THCA; cd ../

