
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





# figure 3 -----

# files from fig3_attributions_gather.R
fig3_data_gathered <- readRDS("bootstrap_analysis/fig3_data_gathered.rds")
fig3_trinuc_data_gathered <- readRDS("bootstrap_analysis/fig3_data_trinuc_gathered.rds")


source("manuscript_analysis_bootstrap/R/fig3_plots_create.R")

fig3_sub_heights <- c(1.5,5)

fig3a1 <- avg_signature_plotter(tumor_type = "LUAD",sig_of_focus = "Tobacco (4,29)",plot_title = "LUAD")
fig3a2 <- variant_attribution_plotter(tumor_type = "LUAD",sig_of_focus = "Tobacco (4,29)")
fig3a <- cowplot::plot_grid(fig3a1,fig3a2,nrow=2,rel_heights = fig3_sub_heights,align = "v")

fig3b1 <- avg_signature_plotter(tumor_type = "LUSC",sig_of_focus = "Tobacco (4,29)",plot_title = "LUSC")
fig3b2 <- variant_attribution_plotter(tumor_type = "LUSC",sig_of_focus = "Tobacco (4,29)")
fig3b <- cowplot::plot_grid(fig3b1,fig3b2,nrow=2,rel_heights = fig3_sub_heights,align = "v")

fig3c1 <- avg_signature_plotter(tumor_type = "SKCMP",sig_of_focus = "UV light (7a–d,38)",plot_title = "Primary SKCM")
fig3c2 <- variant_attribution_plotter(tumor_type = "SKCMP",sig_of_focus = "UV light (7a–d,38)")
fig3c <- cowplot::plot_grid(fig3c1,fig3c2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3d1 <- avg_signature_plotter(tumor_type = "LIHC",sig_of_focus = "Mutagenic chemical exposure (22,24,42,88)",
                                plot_title = "LIHC")
fig3d2 <- variant_attribution_plotter(tumor_type = "LIHC",sig_of_focus = "Mutagenic chemical exposure (22,24,42,88)")
fig3d <- cowplot::plot_grid(fig3d1,fig3d2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3e1 <- avg_signature_plotter(tumor_type = "BLCA",sig_of_focus = "APOBEC (2,13)",plot_title = "BLCA")
fig3e2 <- variant_attribution_plotter(tumor_type = "BLCA",sig_of_focus = "APOBEC (2,13)")
fig3e <- cowplot::plot_grid(fig3e1,fig3e2,nrow=2,rel_heights = fig3_sub_heights,align = "v")

fig3f1 <- avg_signature_plotter(tumor_type = "CESC",sig_of_focus = "APOBEC (2,13)",plot_title = "CESC")
fig3f2 <- variant_attribution_plotter(tumor_type = "CESC",sig_of_focus = "APOBEC (2,13)")
fig3f <- cowplot::plot_grid(fig3f1,fig3f2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3g1 <- avg_signature_plotter(tumor_type = "HNSC_HPVneg",sig_of_focus = "APOBEC (2,13)",plot_title = "HPV-negative HNSC")
fig3g2 <- variant_attribution_plotter(tumor_type = "HNSC_HPVneg",sig_of_focus = "APOBEC (2,13)")
fig3g <- cowplot::plot_grid(fig3g1,fig3g2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3h1 <- avg_signature_plotter(tumor_type = "HNSC_HPVpos",sig_of_focus = "APOBEC (2,13)",plot_title = "HPV-positive HNSC")
fig3h2 <- variant_attribution_plotter(tumor_type = "HNSC_HPVpos",sig_of_focus = "APOBEC (2,13)")
fig3h <- cowplot::plot_grid(fig3h1,fig3h2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3 <- cowplot::plot_grid(fig3a,
                   fig3b,
                   fig3c,
                   fig3d,
                   fig3e,
                   fig3f,
                   fig3g,
                   fig3h,nrow=2,align = "hv",labels = "AUTO")

# cowplot::save_plot(plot = fig3, filename = "manuscript_analysis_bootstrap/figures/fig3.png",base_height = 7.5,base_width = 13)


# +figure3 legend ----- 
test_colors$object <- factor(test_colors$object, levels = test_colors$object)
ggplot(test_colors[4:13,], aes(x=object,y=1,fill=object)) +
  geom_col() +
  scale_fill_manual(values = color_vec,limits=force) +
  theme_classic() + 
  labs(fill="") + 
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.position = "bottom") -> 
  fig3_ggplotlegend

fig3_legend <- suppressWarnings(cowplot::get_legend(fig3_ggplotlegend))




y_label <- ggdraw() +
  draw_label(
    "Variant",
    fontface = 'bold',
    x = 0,
    angle=90,
    y=0.5,
    vjust = 0.5,
    hjust = 0.5,
    size = 15
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

x_label <- ggdraw() +
  draw_label(
    "Proportion of mutational weight or cancer effect attributable to each mutational process",
    fontface = 'bold',
    x = 0.5175,
    angle=0,
    y=0,
    vjust = 1,
    hjust = 0.5,
    size = 15
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(7, 0, 20, 0)
  )



fig3_w_lab <- plot_grid(
  y_label, fig3,
  ncol = 2,
  # rel_heights values control vertical title margins
  rel_widths = c(0.015, 1)
)


fig3_w_lab <- plot_grid(
  fig3_w_lab, x_label,
  nrow = 2,
  # rel_heights values control vertical title margins
  rel_heights = c(1, 0.035)
)


fig3_w_legend <- plot_grid(fig3_legend,fig3_w_lab, ncol=1,rel_heights = c(0.1,1))


cowplot::save_plot(plot = fig3_w_legend, filename = "manuscript_analysis_bootstrap/figures/fig3.png",base_height = 7.5,base_width = 13)





# variant_attribution_plotter(tumor_type = "LUAD",sig_of_focus = "Tobacco (4,29)")
# variant_attribution_plotter(tumor_type = "LUSC",sig_of_focus = "Tobacco (4,29)")
# variant_attribution_plotter(tumor_type = "SKCMP",sig_of_focus = "UV light (7a–d,38)")
# variant_attribution_plotter(tumor_type = "LIHC",sig_of_focus = "Mutagenic chemical exposure (22,24,42,88)")
# variant_attribution_plotter(tumor_type = "BLCA",sig_of_focus = "APOBEC (2,13)")
# variant_attribution_plotter(tumor_type = "CESC",sig_of_focus = "APOBEC (2,13)")
# variant_attribution_plotter(tumor_type = "HNSC_HPVneg",sig_of_focus = "APOBEC (2,13)")
# variant_attribution_plotter(tumor_type = "HNSC_HPVpos",sig_of_focus = "APOBEC (2,13)")
# 
# 
# 
# avg_signature_plotter(tumor_type = "HNSC_HPVpos",sig_of_focus = "APOBEC (2,13)",plot_title = "HNSC HPV negative")
# 























