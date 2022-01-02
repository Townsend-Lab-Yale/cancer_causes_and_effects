
#' Manuscript analyses and figures for Cannataro et al. 
#' Attribution of cancer origins to endogenous, exogenous, and preventable mutational processes


# load in required packages ----- 

library(cancereffectsizeR) # v. 2.3.4
library(tidyverse)
library(data.table)
library(cowplot)
library(RColorBrewer)


# Define manuscript color palette ---- 

source("manuscript_analysis_bootstrap/R/color_palette.R")

# returns color vec
#         color_vec_sbs


# Figure 1, method ---- 
plot_text_size <- 15

# load in specific tumor trinuc data
lusc_trinuc_results <- readRDS(file = "bootstrap_analysis/combined_trinuc_storage_LUSC.rds")


# +Figure 1b, trinuc signature in the tumor ----


source("manuscript_analysis_bootstrap/R/fig1_signatures_in_tumor.R")
fig1_signatures_in_tumor

# returns fig1_signatures_in_tumor, ggplot object 
#         weight_summary, a summary of the bootstrap results of the trinucleotide contexts in the tumor
#         lusc_tumor_trinuc_long, the raw values of the bootstrap results 





# +Figure 1a, trinuc rates within signatures -----

source("manuscript_analysis_bootstrap/R/fig1_signature_definitions.R")
fig1_trinuc_definitions

# returns trinuc_definitions, ggplot object 
#         

# +Figure 1c, probability each source contributed to each variant ----- 

lusc_attribution_results <- readRDS(file = "bootstrap_analysis/combined_attribution_LUSC.rds")


source("manuscript_analysis_bootstrap/R/fig1_sig_contributed_to_variant.R")
fig1_contributions



# +Figure 1d, cancer effect size summary 

lusc_ces <- readRDS("bootstrap_analysis/combined_selection_LUSC.rds")

source("manuscript_analysis_bootstrap/R/fig1_effect_size.R")
fig1_effect_size

# returns fig1_effect_size, ggplot object 


# +Figure 1ef, attributable effect size

source("manuscript_analysis_bootstrap/R/fig1_attributable_effect_size.R")
fig1_attributable_ES
fig1_attributable_effectsize_stacked_nolegend
# fig1_legend


# +Figure1 combine ---- 
fig1a <- cowplot::plot_grid(fig1_trinuc_definitions,labels = "A")
fig1bc <- cowplot::plot_grid(fig1_signatures_in_tumor,fig1_contributions,ncol = 2,labels = c("B","C"))

fig1abc <- cowplot::plot_grid(fig1a,fig1bc,nrow = 2,rel_heights = c(1.2,1))

fig1def <- cowplot::plot_grid(fig1_effect_size,
                              fig1_attributable_ES,
                              fig1_attributable_effectsize_stacked_nolegend,
                              nrow=1,labels = c("D","E","F"),rel_widths = c(1,2,1))

fig1_abcdef <- cowplot::plot_grid(fig1abc,fig1def,nrow=2,rel_heights = c(3,2))
fig1_abcdef_legend <- cowplot::plot_grid(fig1_abcdef,fig1_legend,rel_heights = c(5,1),nrow=2)

cowplot::save_plot(fig1_abcdef_legend,filename = "manuscript_analysis_bootstrap/figures/fig1.png",
                   base_height = 16,base_width = 16)




# figure 2 ----- 




