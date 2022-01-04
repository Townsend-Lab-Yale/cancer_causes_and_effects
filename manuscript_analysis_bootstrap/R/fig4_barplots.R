


fig4_data_df <- rbindlist(fig4_data,idcol = "tumor_type")

fig4_data_df <- rbind(fig4_data_df,fig4_ucec[[1]] %>% mutate(tumor_type ="UCEC"))


fig4_data_df %>% 
  group_by(tumor_type, data_type, signature) %>% 
  summarize(mean_weight = mean(weight)) %>% 
  ungroup() -> 
  mean_signature_weights



mean_signature_weights %>% 
  mutate(signature_process = case_when(
    signature == "SBS1" ~ "Deamination with age, clock-like (1)",
    signature == "SBS5" ~ "Unknown, clock-like (5)",
    signature %in% c("SBS2","SBS13") ~ "APOBEC (2,13)",
    signature %in% c("SBS3") ~ "Defective homologous recombination (3)",
    signature %in% c("SBS4","SBS29") ~ "Tobacco (4,29)",
    signature %in% c("SBS7a","SBS7b","SBS7c","SBS7d","SBS38") ~ "UV light (7a–d,38)",
    signature %in% c("SBS11","SBS31","SBS32","SBS35") ~ "Prior treatment (11,31,32,35)",
    signature %in% c("SBS22","SBS24","SBS42","SBS88") ~ "Mutagenic chemical exposure (22,24,42,88)",
    signature %in% c("SBS16") ~ "Alcohol-associated (16)",
    TRUE ~ "Non-actionable and unknown signatures")) %>% 
  group_by(tumor_type, data_type, signature_process) %>%
  summarize(mean_weight = sum(mean_weight))  %>% 
  ungroup () -> 
  mean_signature_weights

names(color_vec)[6:12]

mean_signature_weights %>% 
  mutate(aging = case_when(
    signature_process %in% names(color_vec)[4:5]~ "clock",
    TRUE ~ "not_clock")) %>% 
  filter(data_type == "effectsize") %>% 
  group_by(tumor_type, aging) %>% 
  summarise(mean_weight = sum(mean_weight)) %>% 
  filter(aging == "clock") %>% 
  arrange(desc(mean_weight)) %>% 
  pull(tumor_type) -> 
  tumor_order


mean_signature_weights$tumor_type <- factor(mean_signature_weights$tumor_type, levels = tumor_order)
mean_signature_weights$signature_process <- factor(mean_signature_weights$signature_process, levels = names(color_vec))

mean_signature_weights$tumor_type <- 
  forcats::fct_recode(mean_signature_weights$tumor_type,
                      `SKCM (Primary)`="SKCMP",
                      `SKCM (Metastases)` = "SKCMM",
                      `HNSC (HPV+)` = "HNSC_HPVpos",
                      `HNSC (HPV−)` = "HNSC_HPVneg",
                      `BRCA (ER–)` = "BRCA_ER_neg",
                      `BRCA (ER+)` = "BRCA_ER_pos")



weight_plot <- ggplot(mean_signature_weights %>% filter(data_type == "trinuc"), 
       aes(x=mean_weight,y=tumor_type,fill=signature_process)) +  
  geom_bar(stat="identity", color="black") + 
  scale_fill_manual(values = color_vec, limits=force) +
  theme_bw() + 
  labs(x="Weight", title = "Mutation") + 
  ylab(NULL) + 
  scale_y_discrete(position = "right") + 
  guides(fill="none") + 
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust=0.5))


axis_text <- ggplot(mean_signature_weights %>% filter(data_type == "trinuc"), 
       aes(x=1,y=tumor_type,fill=signature_process)) +  
  geom_text(aes(label=tumor_type)) +
  theme_bw() + 
  ggtitle("")+ 
  xlab(NULL) +
  # theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  # scale_fill_manual(values = col_vec) + 
  # labs(title="Effects",y="Weight") + 
  # coord_flip() + 
  theme_nothing() +
  # theme(legend.position = "none") + 
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        panel.background=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"))  + 
  # scale_y_reverse() + 
  theme(plot.title = element_text(hjust = 0.5))



effect_plot <- ggplot(mean_signature_weights %>% filter(data_type == "effectsize"), 
                      aes(x=mean_weight,y=tumor_type,fill=signature_process)) +  
  geom_bar(stat="identity", color="black") + 
  scale_fill_manual(values = color_vec, limits=force) +
  theme_bw() + 
  labs(x="Weight", title = "Effects",fill=NULL) + 
  ylab(NULL) + 
  # scale_y_discrete(position = "right") + 
  # guides(fill="none") + 
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  scale_x_reverse()

effects_legend <- suppressWarnings(cowplot::get_legend(plot = effect_plot))






