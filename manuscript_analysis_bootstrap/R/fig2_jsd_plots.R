

# figure 2 bars 

tumors_in_plot <- c("SKCMP",
                    "COAD",
                    "HNSC_HPVneg",
                    "THCA")
# all_median_jsd
# median_quantiles_df
# all_data_for_figure2

median_quantiles_df_forplot <- median_quantiles_df %>%
  filter(tumor_type %in% tumors_in_plot) %>%
  # filter out tumor IDs not plotted 
  filter(!tumor_ID %in% c("G310T",
                          "TCGA-CV-A468-01A-11D-A25Y-08",
                          "TCGA-CV-7568-01A-11D-2229-08",
                          "TCGA-CV-7427-01A-11D-2078-08")) %>%
  # many SKCM have jsd == 0, pick 3 
  mutate(for_plot = case_when(
    tumor_type != "SKCMP" ~ "keep",
    tumor_ID %in% c("TCGA-EB-A82B-01A-11D-A34U-08",
                    "yuchufaT",
                    "TCGA-BF-A1PV-01A-11D-A19A-08",
                    "TCGA-BF-A1PX-01A-12D-A19A-08",
                    "TCGA-BF-A1PZ-01A-11D-A19A-08") ~ "keep",
    TRUE ~ "drop")) %>%
  filter(for_plot == "keep")


all_median_jsd <- all_median_jsd %>%
  mutate(highlighted = case_when(
    tumor_ID %in% median_quantiles_df_forplot$tumor_ID ~ "highlight",
    TRUE ~ "not"))


highlight_colors <- setNames(object = c("white","red"),nm = c("not","highlight"))


jsd_boxplots <- vector(mode = "list",length = length(tumors_in_plot))
names(jsd_boxplots) <- tumors_in_plot



all_median_jsd$tumor_type_label <- forcats::fct_recode(all_median_jsd$tumor_type_label,
                                                       `Primary SKCM`="SKCM, primary",
                                                       `HPV-negative HNSC` = "HNSC (HPV−)")



# build JSD boxplots 
for(tumor_ind in seq_along(jsd_boxplots)){
  
  
  
  this_median_jsd <- all_median_jsd %>%
    filter(tumor_type == tumors_in_plot[tumor_ind]) 
  
  
  
  
  ggplot(data = this_median_jsd, aes(x=median_JSD,y=1)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(data = this_median_jsd %>% 
                  filter(highlighted == "not"), 
                height = 0.1,width =0, aes(fill=highlighted),shape=21) + 
    geom_point(data = this_median_jsd %>% 
                 filter(highlighted == "highlight"),
               aes(fill=highlighted),shape=21) + 
    theme_classic() + 
    scale_fill_manual(values = highlight_colors) + 
    guides(fill = "none") + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + 
    labs(x="Jensen-Shannon Divergence", 
         title= this_median_jsd$tumor_type_label[1]) + 
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = plot_text_size)) -> 
    jsd_boxplots[[tumor_ind]]
  
  
}

# set up the color palette 
all_data_for_figure2 <- all_data_for_figure2 %>%
  mutate(sig_label = case_when(
    signature == "SBS1" ~ "Deamination with age, clock-like (1)",
    signature == "SBS5" ~ "Unknown, clock-like (5)",
    signature %in% c("SBS2_13") ~ "APOBEC (2,13)",
    signature %in% c("SBS3") ~ "Defective homologous recombination (3)",
    signature %in% c("SBS4_29") ~ "Tobacco (4,29)",
    signature %in% c("SBS7_38") ~ "UV light (7a–d,38)",
    signature %in% c("SBStreatment") ~ "Prior treatment (11,31,32,35)",
    signature %in% c("SBSmut") ~ "Mutagenic chemical exposure (22,24,42,88)",
    signature %in% c("SBS16") ~ "Alcohol-associated (16)",
    TRUE ~ "Non-actionable and unknown signatures"))


data_for_figure2 <- all_data_for_figure2 %>%
  filter(tumor_type %in% tumors_in_plot)


data_for_figure2 %>% 
  filter(stringr::str_detect(string = tumor_ID , pattern = "EB-A82B")) %>% 
  arrange(desc(weight)) %>% 
  filter(signature == "SBS1") %>%
  # filter(data_type == "trinuc_weight") %>%
  pull(weight)


data_for_figure2 %>% 
  filter(stringr::str_detect(string = tumor_ID , pattern = "EB-A82B")) %>% 
  arrange(desc(weight)) %>%
  filter(signature %in% c("SBS5","SBS7_38")) %>%
  filter(data_type == "attribution_weight") 
  # pull(weight)


jsd_barplots <- vector(mode = "list", length = length(tumors_in_plot))
names(jsd_barplots) <- tumors_in_plot  

signatures_shown <- NULL

for(tumor_ind in seq_along(tumors_in_plot)){
  
  these_fig2_data <- data_for_figure2 %>%
    filter(tumor_type == tumors_in_plot[tumor_ind]) %>%
    filter(tumor_ID %in% median_quantiles_df_forplot$tumor_ID)
  
  
  tumor_ID_order <- median_quantiles_df_forplot %>%
    filter(tumor_type == tumors_in_plot[tumor_ind]) %>% 
    arrange(desc(quant)) %>% 
    pull(tumor_ID) %>%
    str_sub(1,12)
  
  these_fig2_data$tumor_ID_sub <- str_sub(string =these_fig2_data$tumor_ID, 1,12)
  these_fig2_data$tumor_ID_sub <- factor(these_fig2_data$tumor_ID_sub, 
                                         levels = rev(tumor_ID_order))
  
  
  these_fig2_data <- these_fig2_data %>% 
    filter(weight > 0) %>%
    mutate(x_labels_text = case_when(
      data_type == "trinuc_weight" ~ "SW",
      data_type == "attribution_weight" ~ "CEW")) %>%
  mutate(x_labels_text = factor(x_labels_text,levels=c("SW","CEW")))
  
  
  ggplot(these_fig2_data, 
         aes(x=x_labels_text,y=weight,fill = sig_label)) + 
    geom_bar(stat = "identity",color="black") + 
    theme_classic() + 
    facet_wrap(~tumor_ID_sub,nrow=1) + 
    scale_fill_manual(values = color_vec,limits=force) + 
    guides(fill = "none") + 
    labs(y="Proportion of total weight", x="Type of weight") + 
    theme(text = element_text(size = plot_text_size),
          strip.text.x = element_text(size = plot_text_size-strip_text_smaller)) -> 
    jsd_barplots[[tumor_ind]] 
  
  
  signatures_shown <- rbind(signatures_shown,these_fig2_data)
  
  
  
}



signatures_shown$sig_label <- factor(signatures_shown$sig_label, levels = names(color_vec))

# order for the legend to match figure  

#' 7
#' 1
#' 5
#' unknown
#' 2
#' 3
#' 4
#' 
#' 


color_vec_fig2 <- color_vec[c("UV light (7a–d,38)",
                              "Deamination with age, clock-like (1)",
                              "Unknown, clock-like (5)",
                              "Non-actionable and unknown signatures",
                              "APOBEC (2,13)",
                              "Defective homologous recombination (3)",
                              "Tobacco (4,29)")]
signatures_shown$sig_label <- factor(signatures_shown$sig_label,
                                     levels = names(color_vec_fig2))

signatures_shown %>% 
  ggplot(aes(x=tumor_type,y=weight,fill=sig_label)) + 
  geom_bar(stat="identity", color="black") + 
  scale_fill_manual(values=color_vec_fig2,limits=force) + 
  guides(fill=guide_legend(nrow=3,byrow = T)) +
  theme(legend.position = "bottom") + 
  labs(fill = "") + 
  theme(text = element_text(size = plot_text_size))-> 
  fig2_legend





fig2_legend <- suppressWarnings(cowplot::get_legend(fig2_legend))




