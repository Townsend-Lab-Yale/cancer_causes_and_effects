

# make a big matrix to filter out our tumor of interest
lusc_trinuc_results_rbind <- rbindlist(lusc_trinuc_results,idcol = "bootstrap_sample")

# filter out our tumor of interest
lusc_tumor_trinuc_results <- lusc_trinuc_results_rbind %>%
  filter(Unique_Patient_Identifier == "TCGA-98-A53J-01A-11D-A26M-08")

# make data tidy
lusc_tumor_trinuc_results %>% 
  select(-Unique_Patient_Identifier,-ends_with("snvs"),-group_avg_blended) %>%
  pivot_longer(-bootstrap_sample,names_to = "signature",values_to = "weight") -> 
  lusc_tumor_trinuc_long

# pull out relevant stats
lusc_tumor_trinuc_long %>%
  group_by(signature) %>%
  summarize(qfirst = quantile(weight,0.025),
            median = median(weight),
            mean = mean(weight),
            qlast = quantile(weight, 0.975)) -> 
  weight_summary

# which signatures are in the plot? 
weight_summary %>% 
  filter(qlast > 0) %>% 
  pull(signature) -> 
  signatures_to_plot

# order of the signatures in the figure 
sig_order <- cancereffectsizeR::get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.2")
sig_order <- sig_order$meta$Signature

lusc_tumor_trinuc_long$signature <- factor(lusc_tumor_trinuc_long$signature, levels = sig_order)


weight_summary$signature <- factor(weight_summary$signature,
                                   levels = sig_order)

lusc_tumor_trinuc_long <- lusc_tumor_trinuc_long %>% 
  filter(signature %in% signatures_to_plot)

weight_summary <- weight_summary %>% 
  filter(signature %in% signatures_to_plot)

x_labels <- as.character(sort(weight_summary$signature)) %>%
 stringr::str_remove(pattern = "SBS")

ggplot(weight_summary, aes(x=signature,fill=signature)) + 
  geom_col(aes(y=median),color="black") + 
  geom_jitter(data = lusc_tumor_trinuc_long, aes(y=weight),shape=21,alpha=0.2,width = 0.2,height = 0) +
  geom_errorbar(aes(ymin = qfirst,ymax=qlast), width = 0.3) + 
  theme_classic() + 
  scale_x_discrete(labels = x_labels) + 
  # limits = force gets rid of extra colors
  scale_fill_manual(values = color_vec_sbs, limits=force) + 
  guides(fill = "none") + 
  labs(y="Signature weight",x="Signature") +
  theme(text = element_text(size = plot_text_size)) -> 
  fig1_signatures_in_tumor








