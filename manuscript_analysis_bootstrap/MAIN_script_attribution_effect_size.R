
#' Manuscript analyses and figures for Cannataro et al. 
#' Attribution of cancer origins to endogenous, exogenous, and preventable mutational processes


# load in required packages ----- 

library(cancereffectsizeR) # v. 2.4.0
library(tidyverse)
library(data.table)
library(cowplot)
library(RColorBrewer)
options(dplyr.summarise.inform = FALSE)

# Define manuscript color palette ---- 

source("manuscript_analysis_bootstrap/R/color_palette.R")

# returns color vec
#         color_vec_sbs


# Figure 1, method ---- 
plot_text_size <- 18

# load in specific tumor trinuc data
lusc_trinuc_results <- readRDS(file = "bootstrap_analysis/combined_trinuc_storage_LUSC.rds")


# +Figure 1b, trinuc signature in the tumor ----


source("manuscript_analysis_bootstrap/R/fig1_signatures_in_tumor.R")
# fig1_signatures_in_tumor

# returns fig1_signatures_in_tumor, ggplot object 
#         weight_summary, a summary of the bootstrap results of the trinucleotide contexts in the tumor
#         lusc_tumor_trinuc_long, the raw values of the bootstrap results 



## Called this the first time the script ran, but file is too big for github
# lusc_attribution_results <- readRDS(file = "bootstrap_analysis/combined_attribution_LUSC.rds")
# lusc_attribution_rbind <- rbindlist(lusc_attribution_results,idcol = "bootstrap_sample")
# 
# 
# lusc_attribution_rbind %>%
#   filter(Unique_Patient_Identifier == "TCGA-98-A53J-01A-11D-A26M-08") %>%
#   select(-ends_with("_nrsi"),-Unique_Patient_Identifier) %>%
#   pivot_longer(cols = ends_with("flux"),names_to = "signature",values_to = "weight") ->
#   lusc_attribution_long
# 
# 
# lusc_attribution_rbind %>% 
#   filter(Unique_Patient_Identifier == "TCGA-98-A53J-01A-11D-A26M-08") %>%
#   select(-ends_with("_flux"),-Unique_Patient_Identifier) %>%
#   pivot_longer(cols = ends_with("nrsi"),names_to = "signature",values_to = "weight") %>% 
#   group_by(bootstrap_sample) %>%
#   mutate(weight_prop = weight / sum(weight)) %>% 
#   ungroup() %>%
#   mutate(variant_name = variant_name_extracter(variant)) -> 
#   lusc_attributable_ES_long
# 
# saveRDS(object = lusc_attribution_long, file = "bootstrap_analysis/lusc_attribution_long.rds")
# saveRDS(object = lusc_attributable_ES_long, file = "bootstrap_analysis/lusc_attributable_ES_long.rds")


lusc_attribution_long <- readRDS(file = "bootstrap_analysis/lusc_attribution_long.rds")



# +Figure 1a, trinuc rates within signatures -----

source("manuscript_analysis_bootstrap/R/fig1_signature_definitions.R")
# fig1_trinuc_definitions

# returns trinuc_definitions, ggplot object 
#         

# +Figure 1c, probability each source contributed to each variant ----- 

lusc_attribution_results <- readRDS(file = "bootstrap_analysis/combined_attribution_LUSC.rds")


source("manuscript_analysis_bootstrap/R/fig1_sig_contributed_to_variant.R")
# fig1_contributions



# +Figure 1d, cancer effect size summary 

lusc_ces <- readRDS("bootstrap_analysis/combined_selection_LUSC.rds")

source("manuscript_analysis_bootstrap/R/fig1_effect_size.R")
# fig1_effect_size

# returns fig1_effect_size, ggplot object 


# +Figure 1ef, attributable effect size

lusc_attributable_ES_long <- readRDS(file = "bootstrap_analysis/lusc_attributable_ES_long.rds")

source("manuscript_analysis_bootstrap/R/fig1_attributable_effect_size.R")
# fig1_attributable_ES
# fig1_attributable_effectsize_stacked_nolegend
# fig1_legend


# +Figure1 combine ---- 
fig1a <- cowplot::plot_grid(fig1_trinuc_definitions,labels = "A",label_size = plot_text_size)
fig1bc <- cowplot::plot_grid(fig1_signatures_in_tumor,fig1_contributions,ncol = 2,labels = c("B","C"),label_size = plot_text_size)

fig1abc <- cowplot::plot_grid(fig1a,fig1bc,nrow = 2,rel_heights = c(1.3,1))

fig1def <- cowplot::plot_grid(fig1_effect_size,
                              fig1_attributable_ES,
                              fig1_attributable_effectsize_stacked_nolegend,
                              nrow=1,labels = c("D","E","F"),rel_widths = c(1,2,1),label_size = plot_text_size)

fig1_abcdef <- cowplot::plot_grid(fig1abc,fig1def,nrow=2,rel_heights = c(3.2,2))
fig1_abcdef_legend <- cowplot::plot_grid(fig1_abcdef + theme(plot.margin = unit(c(1,1,-15,1), "mm")),
                                         fig1_legend,rel_heights = c(5,1),nrow=2)

cowplot::save_plot(fig1_abcdef_legend,filename = "manuscript_analysis_bootstrap/figures/fig1.png",
                   base_height = 16,base_width = 16)

cowplot::save_plot(fig1_abcdef_legend,filename = "manuscript_analysis_bootstrap/figures/Cannataro_MBE-21-0913_fig1.eps",
                   base_height = 16,base_width = 16,device=cairo_ps)


# figure 2 ----- 


all_median_jsd <- readRDS(file = "bootstrap_analysis/combined_JSD_medians.rds")

source("manuscript_analysis_bootstrap/R/fig2_median_JSD_analysis.R")
# returns supp_fig2, a ggplot object
#         median_quantiles_df, a dataframe to pull trinuc and attribution data from cluster
#           - contains quantile tumor names and associated data

# supp_fig2



ggplot2::ggsave(filename = "manuscript_analysis_bootstrap/figures/supp_figS2.pdf",plot = supp_fig2,height = 7,width = 7)


# sent median_quantiles_df to the cluster, ran combine_data_for_figure2.R there
# received back all_data_for_figure2.rds



all_data_for_figure2 <- readRDS(file = "bootstrap_analysis/all_data_for_figure2.rds")


strip_text_smaller <- 9
source("manuscript_analysis_bootstrap/R/fig2_jsd_plots.R")
# returns jsd_boxplots
#         jsd_barplots
#         fig2_legend


# +fig2 build -----

heights_dots_v_bars <- c(1.2,3)


fig2a <- cowplot::plot_grid(jsd_boxplots[["SKCMP"]],jsd_barplots[["SKCMP"]],nrow = 2,rel_heights = heights_dots_v_bars)
fig2b <- cowplot::plot_grid(jsd_boxplots[["COAD"]],jsd_barplots[["COAD"]],nrow = 2,rel_heights = heights_dots_v_bars)
fig2c <- cowplot::plot_grid(jsd_boxplots[["HNSC_HPVneg"]],jsd_barplots[["HNSC_HPVneg"]],nrow = 2,rel_heights = heights_dots_v_bars)
fig2d <- cowplot::plot_grid(jsd_boxplots[["THCA"]],jsd_barplots[["THCA"]],nrow = 2,rel_heights = heights_dots_v_bars)

fig2_all <- cowplot::plot_grid(fig2a,fig2b,fig2c,fig2d,nrow = 2,labels = "AUTO",label_size = plot_text_size)

fig2_all_w_legend <- cowplot::plot_grid(fig2_all + theme(plot.margin = unit(c(1,1,-11,1), "mm")),
                                        fig2_legend ,
                                        nrow=2,rel_heights = c(5,1))


cowplot::save_plot(plot = fig2_all_w_legend, 
                   filename = "manuscript_analysis_bootstrap/figures/fig2.png",
                   base_height = 12,
                   base_width = 12)


cowplot::save_plot(plot = fig2_all_w_legend, 
                   filename = "manuscript_analysis_bootstrap/figures/Cannataro_MBE-21-0913_fig2.eps",
                   device = cairo_ps,
                   base_height = 12,
                   base_width = 12)



# figure 3 -----

# files from fig3_attributions_gather.R
fig3_data_gathered <- readRDS("bootstrap_analysis/fig3_data_gathered.rds")
fig3_trinuc_data_gathered <- readRDS("bootstrap_analysis/fig3_data_trinuc_gathered.rds")

plot_text_size <- 20
source("manuscript_analysis_bootstrap/R/fig3_plots_create.R")

fig3_sub_heights <- c(1.6,5)

fig3a1 <- avg_signature_plotter(tumor_type = "LUAD",sig_of_focus = "Tobacco (4,29)",plot_title = "LUAD")
LUAD_var_att <- variant_attribution_plotter(tumor_type = "LUAD",sig_of_focus = "Tobacco (4,29)")
fig3a2 <- LUAD_var_att$this_fig3_attr_plot
fig3a <- cowplot::plot_grid(fig3a1,fig3a2,nrow=2,rel_heights = fig3_sub_heights,align = "v")

fig3b1 <- avg_signature_plotter(tumor_type = "LUSC",sig_of_focus = "Tobacco (4,29)",plot_title = "LUSC")
LUSC_var_att <- variant_attribution_plotter(tumor_type = "LUSC",sig_of_focus = "Tobacco (4,29)")
fig3b2 <- LUSC_var_att$this_fig3_attr_plot
fig3b <- cowplot::plot_grid(fig3b1,fig3b2,nrow=2,rel_heights = fig3_sub_heights,align = "v")

fig3c1 <- avg_signature_plotter(tumor_type = "SKCMP",sig_of_focus = "UV light (7a–d,38)",plot_title = "Primary SKCM")
SKCMP_var_att <- variant_attribution_plotter(tumor_type = "SKCMP",sig_of_focus = "UV light (7a–d,38)")
fig3c2 <- SKCMP_var_att$this_fig3_attr_plot
fig3c <- cowplot::plot_grid(fig3c1,fig3c2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3d1 <- avg_signature_plotter(tumor_type = "LIHC",sig_of_focus = "Mutagenic chemical exposure (22,24,42,88)",
                                plot_title = "LIHC")
LIHC_var_att <- variant_attribution_plotter(tumor_type = "LIHC",sig_of_focus = "Mutagenic chemical exposure (22,24,42,88)")
fig3d2 <- LIHC_var_att$this_fig3_attr_plot
fig3d <- cowplot::plot_grid(fig3d1,fig3d2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3e1 <- avg_signature_plotter(tumor_type = "BLCA",sig_of_focus = "APOBEC (2,13)",plot_title = "BLCA")
BLCA_var_att <- variant_attribution_plotter(tumor_type = "BLCA",sig_of_focus = "APOBEC (2,13)")
fig3e2 <- BLCA_var_att$this_fig3_attr_plot
fig3e <- cowplot::plot_grid(fig3e1,fig3e2,nrow=2,rel_heights = fig3_sub_heights,align = "v")

fig3f1 <- avg_signature_plotter(tumor_type = "CESC",sig_of_focus = "APOBEC (2,13)",plot_title = "CESC")
CESC_var_att <- variant_attribution_plotter(tumor_type = "CESC",sig_of_focus = "APOBEC (2,13)")
fig3f2 <- CESC_var_att$this_fig3_attr_plot
fig3f <- cowplot::plot_grid(fig3f1,fig3f2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3g1 <- avg_signature_plotter(tumor_type = "HNSC_HPVneg",sig_of_focus = "APOBEC (2,13)",plot_title = "HPV-negative HNSC")
HNSC_HPVn_var_att <- variant_attribution_plotter(tumor_type = "HNSC_HPVneg",sig_of_focus = "APOBEC (2,13)")
fig3g2 <- HNSC_HPVn_var_att$this_fig3_attr_plot
fig3g <- cowplot::plot_grid(fig3g1,fig3g2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3h1 <- avg_signature_plotter(tumor_type = "HNSC_HPVpos",sig_of_focus = "APOBEC (2,13)",plot_title = "HPV-positive HNSC")
HNSC_HPVp_var_att <- variant_attribution_plotter(tumor_type = "HNSC_HPVpos",sig_of_focus = "APOBEC (2,13)")
fig3h2 <- HNSC_HPVp_var_att$this_fig3_attr_plot
fig3h <- cowplot::plot_grid(fig3h1,fig3h2,nrow=2,rel_heights = fig3_sub_heights,align = "v")


fig3 <- cowplot::plot_grid(fig3a,
                   fig3b,
                   fig3c,
                   fig3d,
                   fig3e,
                   fig3f,
                   fig3g,
                   fig3h,nrow=4,align = "hv",labels = "AUTO",label_size = plot_text_size)

# cowplot::save_plot(plot = fig3, filename = "manuscript_analysis_bootstrap/figures/fig3.png",base_height = 7.5,base_width = 13)


# +figure3 legend ----- 
test_colors$object <- factor(test_colors$object, levels = c(
  "Tobacco (4,29)",
  "Non-actionable and unknown signatures",
  "Defective homologous recombination (3)",
  "Unknown, clock-like (5)",
  "Deamination with age, clock-like (1)",
  "UV light (7a–d,38)",
  "APOBEC (2,13)",
  "Mutagenic chemical exposure (22,24,42,88)",
  "Alcohol-associated (16)")
  )

test_colors <- test_colors %>%
  filter(!is.na(object))

color_vec_fig3 <- color_vec[levels(test_colors$object)]

ggplot(test_colors, aes(x=object,y=1,fill=object)) +
  geom_col() +
  scale_fill_manual(values = color_vec_fig3,limits=force) +
  theme_classic() + 
  labs(fill="") + 
  guides(fill=guide_legend(nrow=5,byrow = T)) +
  theme(legend.position = "bottom") + 
  theme(text = element_text(size = plot_text_size),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0, 'cm'))-> 
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
    size = plot_text_size
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 15, 0, 15)
  )

x_label <- ggdraw() +
  draw_label(
    "Proportion of mutational weight or cancer effect \nattributable to each mutational process",
    fontface = 'bold',
    x = 0.5175,
    angle=0,
    y=0,
    vjust = 1,
    hjust = 0.5,
    size = plot_text_size
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(7, 0, 50, 0)
  )



fig3_w_lab <- plot_grid(
  y_label, fig3,
  ncol = 2,
  # rel_heights values control vertical title margins
  rel_widths = c(0.020, 1)
)


fig3_w_lab <- plot_grid(
  fig3_w_lab, x_label,
  nrow = 2,
  # rel_heights values control vertical title margins
  rel_heights = c(1, 0.035)
)


fig3_w_legend <- plot_grid(fig3_legend,fig3_w_lab, ncol=1,rel_heights = c(0.1,1))


cowplot::save_plot(plot = fig3_w_legend, filename = "manuscript_analysis_bootstrap/figures/fig3.png",base_height = 20,base_width = 10)

cowplot::save_plot(plot = fig3_w_legend, filename = "manuscript_analysis_bootstrap/figures/Cannataro_MBE-21-0913_fig3.eps",
                   base_height = 20,base_width = 10,device= cairo_ps)


# +supp table relating to figure 3 data ------ 

fig3_supp_table <- rbind(LUAD_var_att$fig3_quantile_data_all,
      LUSC_var_att$fig3_quantile_data_all,
      SKCMP_var_att$fig3_quantile_data_all,
      LIHC_var_att$fig3_quantile_data_all,
      BLCA_var_att$fig3_quantile_data_all,
      CESC_var_att$fig3_quantile_data_all,
      HNSC_HPVn_var_att$fig3_quantile_data_all,
      HNSC_HPVp_var_att$fig3_quantile_data_all)

write_csv(x = fig3_supp_table,file = "manuscript_analysis_bootstrap/tables/table_S1_fig3supp_info.csv")

# figure 4 ------ 


fig4_data <- readRDS(file = "bootstrap_analysis/fig_4_storage.rds")
fig4_ucec <- readRDS(file = "bootstrap_analysis/fig_4_storage_UCEC.rds")

source("manuscript_analysis_bootstrap/R/fig4_dotplot.R")


write_csv(x = weight_and_effect_data_forSuppTable,file = "manuscript_analysis_bootstrap/tables/table_S2_weight_compare.csv")



# fig4dotplot

source("manuscript_analysis_bootstrap/R/fig4_barplots.R")

fig4_bars <- cowplot::plot_grid(weight_plot, axis_text, effect_plot + guides(fill="none"),
                                nrow = 1,rel_widths = c(1,.6,1,.5),align = "h",labels = c("B","","C"))

fig4_bars_w_legend <- cowplot::plot_grid(fig4_bars,effects_legend,nrow = ,rel_widths = c(1,.5))

fig4 <- cowplot::plot_grid(fig4dotplot,fig4_bars_w_legend,nrow = 2,rel_heights = c(2,1),labels = c("A",""))

cowplot::save_plot(plot = fig4,filename = "manuscript_analysis_bootstrap/figures/fig4.png",base_height = 12,base_width = 9)

cowplot::save_plot(plot = fig4,filename = "manuscript_analysis_bootstrap/figures/Cannataro_MBE-21-0913_fig4.eps",
                   base_height = 12,base_width = 9,device=cairo_ps)

# supp figure 1 ------

# file created with variant_number_extract.R on the cluster
variant_numbers <- readRDS("bootstrap_analysis/variant_numbers.rds")


source("manuscript_analysis_bootstrap/R/suppfig1_variant_numbers.R")
# proportion_of_tumors_greater_than_50_subs

ggplot2::ggsave(filename = "manuscript_analysis_bootstrap/figures/figure_s1_snv_50ormore.pdf",
                plot= proportion_of_tumors_greater_than_50_subs,height = 4,width = 6)


# supp table YG tumors removed ---- 
HNSC_HPVneg_rt <- get(load("input_data/from_cluster/HNSC_HPVneg_removed_tumors.RData"))
SKCM_prim_rt <- get(load("input_data/from_cluster/SKCMP_removed_tumors.RData"))
SKCM_met_rt <- get(load("input_data/from_cluster/SKCMM_removed_tumors.RData"))

removed_yg_tumors <- unique(c(HNSC_HPVneg_rt,SKCM_prim_rt,SKCM_met_rt))


write.table(x = removed_yg_tumors,file = "manuscript_analysis_bootstrap/tables/table_S3_removed_YG_tumors.csv",quote = F,row.names = F,col.names = F)


# in-text descriptions ----- 

weight_and_effect_data %>% 
  ungroup() %>%
  # filter(signature %in% c("SBS2","SBS13")) %>% 
  pivot_wider(values_from = mean_weight,names_from = data_type) %>% 
  mutate(effect_higher = case_when(
    effectsize > trinuc ~ TRUE, 
    TRUE ~ FALSE)) -> 
  compare_effect_v_trinuc

message("APOBEC contributes more to mutation than effect size: ")
compare_effect_v_trinuc %>% 
  filter(signature %in% c("SBS2","SBS13")) %>% 
  count(significant, effect_higher) %>% 
  print()

message("Signatures where mutational weights are consistently and significantly higher than 
        cancer effects:")
compare_effect_v_trinuc %>% 
  group_by(signature) %>% 
  count(significant,effect_higher) %>% 
  arrange(desc(n)) %>% 
  filter(effect_higher == F,significant == T) %>%
  print(n=15)

message("Signatures where cancer effects are consistently and significantly higher than 
        trinucs:")
compare_effect_v_trinuc %>% 
  group_by(signature) %>% 
  count(significant,effect_higher) %>% 
  arrange(desc(n)) %>% 
  filter(effect_higher == T,significant == T) %>%
  print(n=15)



message("Age associated signature 1: ")
compare_effect_v_trinuc %>% 
  filter(signature %in% c("SBS1")) %>% 
  count(significant, effect_higher) %>% 
  print()


message("Age associated signature 5: ")
compare_effect_v_trinuc %>% 
  filter(signature %in% c("SBS5")) %>% 
  count(significant, effect_higher) %>% 
  print()


message("Aging in LGG:")
compare_effect_v_trinuc %>% 
  filter(tumor_type == "LGG") %>%
  filter(signature %in% c("SBS1","SBS5")) %>% 
  print()


message("APOBEC in THCA:")
compare_effect_v_trinuc %>% 
  filter(tumor_type == "THCA") %>%
  filter(signature %in% c("SBS2","SBS13")) %>% 
  print()


message("Aging in THCA:")
compare_effect_v_trinuc %>% 
  filter(tumor_type == "THCA") %>%
  filter(signature %in% c("SBS1","SBS5")) %>% 
  print()



message("APOBEC in THCA mutation weight")
mean_signature_weights %>% 
  filter(tumor_type == "THCA") %>%
  filter(data_type == "trinuc") %>% 
  print(n=10)


message("APOBEC in THCA cancer effect weight")
mean_signature_weights %>% 
  filter(tumor_type == "THCA") %>%
  filter(data_type == "effectsize") %>% 
  print(n=10)

fig4_stats %>%
  filter(tumor_type == "THCA") %>%
  filter(signature %in% c("SBS2","SBS13")) %>%
  print()


message("Mutational effect vs cancer effect tobacco in lung cancer:")
mean_signature_weights %>% 
  filter(tumor_type %in% c("LUAD","LUSC")) %>%
  filter(signature_process == "Tobacco (4,29)") %>% 
  print(n=10)


compare_effect_v_trinuc %>% 
  filter(tumor_type %in% c("LUAD","LUSC")) %>%
  filter(signature %in% c("SBS4","SBS29")) %>% 
  print(n=10)

fig4_stats %>%
  filter(tumor_type %in% c("LUAD","LUSC")) %>%
  filter(signature %in% c("SBS4","SBS29")) %>%
  print()



message("Signature 1 among cancers:")
mean_signature_weights %>%
  filter(signature_process == "Deamination with age, clock-like (1)") %>% 
  filter(data_type == "effectsize") %>% 
  arrange(mean_weight) %>% 
  print(n=Inf)



message("Signature 1 and 5 among cancers:")
mean_signature_weights %>%
  filter(signature_process %in% c("Deamination with age, clock-like (1)","Unknown, clock-like (5)")) %>% 
  filter(data_type == "effectsize") %>% 
  group_by(tumor_type) %>% 
  summarize(total_mean_weight = sum(mean_weight)) %>%
  arrange(total_mean_weight) -> all_age_associated

all_age_associated %>%
  print(n=Inf)



weight_and_effect_data

sig_context <- cancereffectsizeR::get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.2")
etiol <- sig_context$meta$Etiology
names(etiol) <- sig_context$meta$Signature
  
weight_and_effect_data$etiol <- etiol[as.character(weight_and_effect_data$signature)]

message("Unknown signatures cancer effect size among cancers:")
weight_and_effect_data %>%
  filter(str_detect(string = etiol,pattern = "Unknown")) %>%
  filter(!etiol == "Unknown, clock-like") %>%
  filter(!signature == "SBS16") %>%
  filter(data_type == "effectsize") %>%
  group_by(tumor_type) %>% 
  summarize(total_unknown = sum(mean_weight)) %>%
  arrange(desc(total_unknown))
