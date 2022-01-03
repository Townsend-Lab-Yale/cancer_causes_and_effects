# load in JSD files from cluster 

# build with combining_jsd_median_files.R 

# all_median_jsd <- readRDS(file = "bootstrap_analysis/combined_JSD_medians.rds")


# all_median_jsd %>% count(tumor_type) 

all_median_jsd %>%
  group_by(tumor_type) %>%
  summarize(group_median = median(median_JSD)) %>%
  arrange(desc(group_median)) %>% 
  pull(tumor_type) -> 
  tumor_order
  

all_median_jsd$tumor_type <- factor(all_median_jsd$tumor_type, levels = tumor_order)

all_median_jsd$tumor_type_label <- forcats::fct_recode(all_median_jsd$tumor_type,
                                                              `SKCM, primary`="SKCMP",
                                                              `SKCM, metastases` = "SKCMM",
                                                              `HNSC (HPV+)` = "HNSC_HPVpos",
                                                              `HNSC (HPV−)` = "HNSC_HPVneg",
                                                              `BRCA (ER–)` = "BRCA_ER_neg",
                                                              `BRCA (ER+)` = "BRCA_ER_pos")


ggplot(all_median_jsd, aes(x=tumor_type_label,y=median_JSD)) + 
  geom_jitter(width=0.25,alpha=0.5) + 
  geom_boxplot(outlier.shape = NA,alpha=0.6) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + 
  labs(x="Tumor type", y="Jensen-Shannon Divergence") -> 
  supp_fig2




# which tumors are in the quantiles we want? 

median_jsd_list <- split(x = all_median_jsd, f = all_median_jsd$tumor_type)

# median_tumor_type <- median_jsd_list$HNSC_HPVneg # debug mode

find_quantiles <- function(median_tumor_type){
  
  median_tumor_type <- median_tumor_type %>% arrange(desc(median_JSD))
  
  values_on_quantiles <- median_tumor_type$median_JSD[c(1,
                                              round(nrow(median_tumor_type)/4),
                                              round(nrow(median_tumor_type)/2),
                                              round(nrow(median_tumor_type)/2) + round(nrow(median_tumor_type)/4),
                                              nrow(median_tumor_type))]
  
  values_on_quant_tib <- tibble(values = values_on_quantiles, quant = c(1,0.75,0.5,0.25,0))
  
  
  median_quantiles <- median_tumor_type[which(median_tumor_type$median_JSD %in% values_on_quantiles),]
  
  median_quantiles$quant <- as.numeric(NA)
  for(row_ind in 1:nrow(values_on_quant_tib)){
    
    median_quantiles[which(median_quantiles$median_JSD == pull(values_on_quant_tib[row_ind,"values"])),"quant"] <- 
      values_on_quant_tib[row_ind,"quant"]
    
    
  }
  
  
  return(median_quantiles)
  
}


median_quantiles_df <- purrr::map_dfr(median_jsd_list,find_quantiles)


saveRDS(object = median_quantiles_df, file = "scripts_for_cluster/with_bootstrap/tumors_for_JSD_plot.rds")




