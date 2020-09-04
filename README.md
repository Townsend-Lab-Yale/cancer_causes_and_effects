A repository associated with the analyses presented within Cannataro et al. Attributing cancer 
to endogenous, exogenous, and actionable mutational processes. 

Cancer effect sizes are calculated as described within https://townsend-lab-yale.github.io/cancereffectsizeR/ . 

The RMarkdown script `Population_scaled_analysis_manuscript_calculations.Rmd` contains the code to calculate our main findings and generate all of the figures in our manuscript. 

Analyses that generate input to `Population_scaled_analysis_manuscript_calculations.Rmd` were run on the high-performance computing cluster, for specific files and parameters used, see `scripts_for_cluster/selection_and_pop_level_cluster_submitter_batch.sh`. This shell script ultimately calls `scripts_for_cluster/selection_and_pop_cluster.R`, which performs the effect size, mutational source effect size, and Jensen-Shannon Divergence calculations. 










