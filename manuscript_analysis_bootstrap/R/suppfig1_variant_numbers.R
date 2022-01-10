

# supp fig 1

# # file created with variant_number_extract.R on the cluster
# variant_numbers <- readRDS("bootstrap_analysis/variant_numbers.rds")



less_than_50 <- variant_numbers %>% 
  group_by(tumor_name) %>% 
  summarize(less_than_50 = length(which(n < 50))/ n()) %>% 
  arrange((less_than_50))



less_than_50 <- less_than_50 %>%
  mutate(tumor_type = tumor_name) 
  

less_than_50$tumor_type <- factor(less_than_50$tumor_type, levels=less_than_50$tumor_type)
less_than_50$tumor_type <- forcats::fct_recode(less_than_50$tumor_type,
                                               `SKCM, primary`="SKCM_primary",
                                               `SKCM, metastases` = "SKCM_metastasis",
                                               `HNSC (HPV+)` = "HNSC_HPVpos",
                                               `HNSC (HPV−)` = "HNSC_HPVneg",
                                               `BRCA (ER–)` = "BRCA_ER_neg",
                                               `BRCA (ER+)` = "BRCA_ER_pos")

total_tumors <- variant_numbers %>% 
  group_by(tumor_name) %>%
  summarize(total_tumors_n = n_distinct(Unique_Patient_Identifier))

total_tumors <- left_join(x = total_tumors,y=less_than_50,by=c("tumor_name"))



proportion_of_tumors_greater_than_50_subs <- ggplot(less_than_50) +
  geom_point(aes(x=tumor_type,y=1-less_than_50)) + 
  geom_text(data=total_tumors, aes(x=tumor_type,y=(1-less_than_50)+.03,label=total_tumors_n),angle=90,hjust=0,size=3.5) +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) + 
  coord_cartesian(ylim=c(0,1.15)) + 
  scale_y_continuous(breaks=c(0,.2,.4,.6,.8,1)) + 
  labs(x="Tumor type", y="Proportion of tumors with \ngreater than 50 substitutions")







