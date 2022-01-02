# load in JSD files from cluster 

# build with combining_jsd_median_files.R 

all_median_jsd <- readRDS(file = "bootstrap_analysis/combined_JSD_medians.rds")


# all_median_jsd %>% count(tumor_type) 

all_median_jsd %>%
  group_by(tumor_type) %>%
  summarize(group_median = median(median_JSD)) %>%
  arrange(desc(group_median)) %>% 
  pull(tumor_type) -> 
  tumor_order
  

all_median_jsd$tumor_type <- factor(all_median_jsd$tumor_type, levels = tumor_order)

all_median_jsd$tumor_type <- forcats::fct_recode(all_median_jsd$tumor_type,
                                                              `SKCM, primary`="SKCM_primary",
                                                              `SKCM, metastases` = "SKCM_metastasis",
                                                              `HNSC (HPV+)` = "HNSC_HPVpos",
                                                              `HNSC (HPV−)` = "HNSC_HPVneg",
                                                              `BRCA (ER–)` = "BRCA_ER_neg",
                                                              `BRCA (ER+)` = "BRCA_ER_pos")


ggplot(all_median_jsd, aes(x=tumor_type,y=median_JSD)) + 
  geom_boxplot() +
  geom_jitter(width=0.25,alpha=0.5) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + 
  labs(x="Tumor type", y="Jensen-Shannon Divergence")










