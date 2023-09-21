A repository associated with [Attribution of Cancer Origins to Endogenous, Exogenous, and Preventable Mutational Processes](https://academic.oup.com/mbe/article/39/5/msac084/6570859), *Molecular and Biological Evolution* (Cannataro, Mandell, and Townsend 2023).

**New:** The [cancereffectsizeR (https://townsend-lab-yale.github.io/cancereffectsizeR/)] package now has a [built-in function ](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/mutational_signature_effects.html) to attribute cancer effects to mutational processes using our method. Naturally, for your own analyses, you should use this actively maintained package rather than extracting code from this repository. :)



The latest iteration of our manuscript uses a bootstrap approach to resample variants within tumors, and subsequently generate confidence intervals around detected signatures and downstream analyses.  The manuscript figures and analyses are generated within [`manuscript_analysis_bootstrap/MAIN_script_attribution_effect_size.R`](https://github.com/Townsend-Lab-Yale/cancer_causes_and_effects/blob/master/manuscript_analysis_bootstrap/MAIN_script_attribution_effect_size.R). This script sources various R scripts in the [`R`](https://github.com/Townsend-Lab-Yale/cancer_causes_and_effects/tree/master/manuscript_analysis_bootstrap/R) subdirectory to generate perform calculations and generate figures. It also calls input data from `bootstrap_analysis/`, which were generated on the Yale HPC clusters. The scripts to generate these data, gather them, and retrieve them are stored within [`scripts_for_cluster/with_bootstrap`](https://github.com/Townsend-Lab-Yale/cancer_causes_and_effects/tree/master/scripts_for_cluster/with_bootstrap) 

Before bootstrap approach: 

The RMarkdown script `Population_scaled_analysis_manuscript_calculations.Rmd` contains the code to calculate our main findings and generate all of the figures in our manuscript. 

Analyses that generate input to `Population_scaled_analysis_manuscript_calculations.Rmd` were run on the high-performance computing cluster, for specific files and parameters used, see `scripts_for_cluster/selection_and_pop_level_cluster_submitter_batch.sh`. This shell script ultimately calls `scripts_for_cluster/selection_and_pop_cluster.R`, which performs the effect size, mutational source effect size, and Jensen-Shannon Divergence calculations. 










