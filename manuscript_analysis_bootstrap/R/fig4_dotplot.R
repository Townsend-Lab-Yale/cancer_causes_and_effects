

# figure 4 data combine

fig4_data_df <- rbindlist(fig4_data,idcol = "tumor_type")

fig4_data_df <- rbind(fig4_data_df,fig4_ucec[[1]] %>% mutate(tumor_type ="UCEC"))
# 
# fig4_data_df %>% 
#   group_by(tumor_type,signature) %>% 
#   summarize(total_weight = sum(weight)) %>% 
#   filter(total_weight == 0) -> 
#   to_drop
# 
# fig4_data_df <- anti_join(x = fig4_data_df,y = to_drop,by = c("tumor_type","signature"))

# just samples with 30 or more >0 incidence per tumor type for stats

fig4_data_df %>% 
  filter(data_type == "trinuc") %>% 
  filter(weight > 0) %>%
  # group_by(tumor_type, signature) %>%
  count(tumor_type, signature) %>%
  filter(n<30) -> 
  not_enough_samples
  
fig4_data_df <- anti_join(x = fig4_data_df, y = not_enough_samples, by=c("tumor_type","signature"))

fig4_data_df %>%
  filter((data_type == "trinuc" & weight == 0)) -> 
  no_weight

fig4_data_df <- anti_join(x = fig4_data_df, y = no_weight, by=c("Unique_Patient_Identifier","signature"))

fig4_stats <- fig4_data_df %>% 
  # filter(!Unique_Patient_Identifier %in% no_weight) %>%
  group_by(tumor_type, signature) %>% 
  summarize(results = list(broom::tidy(wilcox.test(weight ~ data_type, paired=T)))) %>% 
  unnest(cols = c(results)) %>% 
  mutate(significant = ifelse(p.value < 0.05,T,F))



weight_and_effect_data <- fig4_data_df %>% 
  group_by(tumor_type, signature, data_type) %>%
  summarize(mean_weight = mean(weight), 
            lower_conf = Rmisc::CI(weight)['lower'], 
            upper_conf = Rmisc::CI(weight)['upper'])


weight_and_effect_data_forSuppTable <- weight_and_effect_data
  
weight_and_effect_data$significant <- NA


for(row_ind in 1:nrow(fig4_stats)){
  
  weight_and_effect_data[
    which(weight_and_effect_data$tumor_type == fig4_stats$tumor_type[row_ind] & 
          weight_and_effect_data$signature == fig4_stats$signature[row_ind]),"significant"
  ] <- fig4_stats$significant[row_ind]
  
}

weight_and_effect_data$tumor_type <- 
  forcats::fct_recode(weight_and_effect_data$tumor_type,
                      `SKCM (Primary)`="SKCMP",
                      `SKCM (Metastases)` = "SKCMM",
                      `HNSC (HPV+)` = "HNSC_HPVpos",
                      `HNSC (HPV−)` = "HNSC_HPVneg",
                      `BRCA (ER–)` = "BRCA_ER_neg",
                      `BRCA (ER+)` = "BRCA_ER_pos")

v_lines <- seq(1,length(unique(weight_and_effect_data$tumor_type)),1) + 0.5
h_lines <- seq(1,length(unique(weight_and_effect_data$signature)),1) + 0.5

sig_levels <- cancereffectsizeR::get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.2")
sig_levels <- sig_levels$meta$Signature

weight_and_effect_data$signature <- factor(weight_and_effect_data$signature,levels = rev(sig_levels))

sigs_present <- as.character(unique(weight_and_effect_data$signature))

sigs_present <- factor(sigs_present, levels=sig_levels[sig_levels %in% sigs_present])


fig4dotplot <- ggplot() + 
  geom_text(data = subset(weight_and_effect_data,data_type=="trinuc" & significant == FALSE), 
            aes(x=tumor_type, y=signature,size=mean_weight),
            label="\u25D6", family = "Arial Unicode MS",col="gray50") + 
  geom_text(data = subset(weight_and_effect_data,data_type=="effectsize" & significant == FALSE), 
            aes(x=tumor_type, y=signature,size=mean_weight),
            label="\u25D7", family = "Arial Unicode MS",col="gray90") + 
  geom_text(data = subset(weight_and_effect_data,data_type=="trinuc" & significant == TRUE), 
            aes(x=tumor_type, y=signature,size=mean_weight),
            label="\u25D6", family = "Arial Unicode MS") + 
  geom_text(data = subset(weight_and_effect_data,data_type=="effectsize" & significant == TRUE), 
            aes(x=tumor_type, y=signature,size=mean_weight),
            label="\u25D7", family = "Arial Unicode MS",col="red") + 
  theme_classic()  + 
  scale_x_discrete(position = "top") + 
  scale_y_discrete(limits=rev(levels(sigs_present)),
    labels=gsub(pattern = "SBS",replacement = "",x = rev(levels(sigs_present))))  +
  theme(axis.text.x = element_text(angle=55,hjust=0),
        text=element_text(size=plot_text_size-5)) + 
  labs(x="Tumor type",
       y="Signature") + 
  scale_size(range = c(0, 10),breaks = c(0,0.05,0.1,0.2,0.4,0.6,1),name = "Proportional signature weight | Cancer effect weight") +
  guides(size=guide_legend(nrow=3,byrow=TRUE)) + 
  geom_vline(xintercept = v_lines,alpha=0.1) +
  theme(legend.margin = margin(-0.5,unit="cm")) + 
  geom_hline(yintercept = h_lines,alpha=0.1) +
  theme(legend.position = "bottom")  







# weight_and_effect_data


# fig4_data_df %>% 
#   group_by(tumor_type, signature,data_type) %>%
#   summarize(mean_weight = mean(weight), 
#             lower_conf = Rmisc::CI(weight)['lower'], 
#             upper_conf = Rmisc::CI(weight)['upper'])
  








weight_and_effect_data_forSuppTable <- weight_and_effect_data_forSuppTable %>% 
  ungroup() %>%
  mutate(data_type = case_when(
    data_type == "trinuc" ~ "Proportional signature weight", 
    data_type == "effectsize" ~ "Cancer effect weight")) %>%
  # pivot_wider(values_from = c("mean_weight","lower_conf","upper_conf"),names_from = data_type) %>%
  # dplyr::rename(
  #   mean_effectsize_weight = effectsize,
  #   mean_trinuc_weight = trinuc) %>% 
  left_join(fig4_stats, by=c("tumor_type","signature"))


