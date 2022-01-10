A repository associated with the analyses presented within Cannataro et al. Attributing cancer 
to endogenous, exogenous, and actionable mutational processes. 

Cancer effect sizes are calculated as described within https://townsend-lab-yale.github.io/cancereffectsizeR/ . 

After calculating effect size using the pipeline outlined at the above link, you can feed the cancereffectsizeR object output into the [`population_scaled_effect_per_tumor`](https://github.com/Townsend-Lab-Yale/cancer_causes_and_effects/blob/master/R/population_scaled_effect_per_tumor.R) function to calculate both the probability a variant was the result of a specific mutational signature (columns ending in _flux) or the attributed effect size (columns ending in _nrsi). 

The latest iteration of our manuscript (in review) uses bootstrap approach to resample variants within tumors, and subsequently generate confidence intervals around detected signatures and downstream analyses.  The manuscript figures and analyses are generated within [`manuscript_analysis_bootstrap/MAIN_script_attribution_effect_size.R`](https://github.com/Townsend-Lab-Yale/cancer_causes_and_effects/blob/master/manuscript_analysis_bootstrap/MAIN_script_attribution_effect_size.R). This script sources various R scripts in the [`R`](https://github.com/Townsend-Lab-Yale/cancer_causes_and_effects/tree/master/manuscript_analysis_bootstrap/R) subdirectory to generate perform calculations and generate figures. It also calls input data from `bootstrap_analysis/`, which were generated on the Yale HPC clusters. The scripts to generate these data, gather them, and retrieve them are stored within [`scripts_for_cluster/with_bootstrap`](https://github.com/Townsend-Lab-Yale/cancer_causes_and_effects/tree/master/scripts_for_cluster/with_bootstrap) 

Before bootstrap approach: 

The RMarkdown script `Population_scaled_analysis_manuscript_calculations.Rmd` contains the code to calculate our main findings and generate all of the figures in our manuscript. 

Analyses that generate input to `Population_scaled_analysis_manuscript_calculations.Rmd` were run on the high-performance computing cluster, for specific files and parameters used, see `scripts_for_cluster/selection_and_pop_level_cluster_submitter_batch.sh`. This shell script ultimately calls `scripts_for_cluster/selection_and_pop_cluster.R`, which performs the effect size, mutational source effect size, and Jensen-Shannon Divergence calculations. 










