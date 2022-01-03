
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


all_median_jsd <- readRDS(file = "bootstrap_analysis/combined_JSD_medians.rds")

source("manuscript_analysis_bootstrap/R/fig2_median_JSD_analysis.R")
# returns supp_fig2, a ggplot object
#         median_quantiles_df, a dataframe to pull trinuc and attribution data from cluster
#           - contains quantile tumor names and associated data

supp_fig2



ggplot2::ggsave(filename = "manuscript_analysis_bootstrap/figures/supp_figS2.png",plot = supp_fig2,height = 7,width = 7)


# sent median_quantiles_df to the cluster, ran combine_data_for_figure2.R there
# received back all_data_for_figure2.rds



all_data_for_figure2 <- readRDS(file = "bootstrap_analysis/all_data_for_figure2.rds")

source("manuscript_analysis_bootstrap/R/fig2_jsd_plots.R")
# returns jsd_boxplots
#         jsd_barplots
#         fig2_legend


# +fig2 build -----

fig2a <- cowplot::plot_grid(jsd_boxplots[["SKCMP"]],jsd_barplots[["SKCMP"]],nrow = 2,rel_heights = c(1,3))
fig2b <- cowplot::plot_grid(jsd_boxplots[["COAD"]],jsd_barplots[["COAD"]],nrow = 2,rel_heights = c(1,3))
fig2c <- cowplot::plot_grid(jsd_boxplots[["HNSC_HPVneg"]],jsd_barplots[["HNSC_HPVneg"]],nrow = 2,rel_heights = c(1,3))
fig2d <- cowplot::plot_grid(jsd_boxplots[["THCA"]],jsd_barplots[["THCA"]],nrow = 2,rel_heights = c(1,3))

fig2_all <- cowplot::plot_grid(fig2a,fig2b,fig2c,fig2d,nrow = 2,labels = "AUTO")

fig2_all_w_legend <- cowplot::plot_grid(fig2_all,fig2_legend,nrow=2,rel_heights = c(5,1))

cowplot::save_plot(plot = fig2_all_w_legend, 
                   filename = "manuscript_analysis_bootstrap/figures/fig2.png",
                   base_height = 12,
                   base_width = 12)



















