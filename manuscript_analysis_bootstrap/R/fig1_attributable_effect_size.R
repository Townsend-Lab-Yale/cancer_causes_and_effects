
# attributable effect size

lusc_attribution_rbind %>% #from fig1_sig_contributed_to_variant.R
  filter(Unique_Patient_Identifier == "TCGA-98-A53J-01A-11D-A26M-08") %>%
  select(-ends_with("_flux"),-Unique_Patient_Identifier) %>%
  pivot_longer(cols = ends_with("nrsi"),names_to = "signature",values_to = "weight") %>% 
  group_by(bootstrap_sample) %>%
  mutate(weight_prop = weight / sum(weight)) %>% 
  ungroup() %>%
  mutate(variant_name = variant_name_extracter(variant)) -> 
  lusc_attributable_ES_long


lusc_attributable_ES_long$variant_name <- factor(lusc_attributable_ES_long$variant_name, levels = c("TP53 R282W",
                                                                                            "KLF5 E419Q",
                                                                                            "OR2T34 L163L"))

lusc_attributable_ES_long <- lusc_attributable_ES_long %>%
  mutate(sig_label = stringr::str_remove(string = signature, pattern = "_nrsi"))

lusc_attributable_ES_long$sig_label <- factor(lusc_attributable_ES_long$sig_label, levels = sig_order)


lusc_attributable_ES_long %>%
  group_by(sig_label,variant_name) %>%
  summarize(qfirst = quantile(weight_prop,0.025),
            median = median(weight_prop),
            mean = mean(weight_prop),
            qlast = quantile(weight_prop, 0.975)) -> 
  attribution_summary

attribution_summary %>%
  filter(qlast>0) %>%
  pull(sig_label) %>%
  as.character() -> 
  sigs_to_plot

lusc_attributable_ES_long <- lusc_attributable_ES_long %>%
  filter(sig_label %in% sigs_to_plot)


# lusc_attribution_long %>% 
#   ggplot(aes(x=sig_label,y=weight,fill = sig_label)) + 
#   geom_jitter(shape = 21, alpha=0.2) +
#   scale_fill_manual(values  = color_vec_sbs,limits = force) + 
#   facet_wrap(~variant_name) + 
#   theme_classic() + 
#   theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) + 
#   guides(fill="none")

attribution_summary <- attribution_summary %>%
  filter(sig_label %in%  sigs_to_plot)

x_labels <- as.character(sort(unique(attribution_summary$sig_label))) %>%
  stringr::str_remove(pattern = "SBS")

borders_colors <- color_vec[levels(attribution_summary$variant_name)]


ggplot(attribution_summary, aes(x=sig_label,fill=sig_label)) + 
  geom_col(aes(y=median),color="black") + 
  geom_jitter(data = lusc_attributable_ES_long, aes(y=weight_prop), shape = 21, alpha = 0.2, width = 0.2,height = 0) + 
  geom_errorbar(aes(ymin = qfirst, ymax = qlast), width = 0.3) + 
  theme_classic() + 
  scale_x_discrete(labels = x_labels) + 
  facet_wrap(~variant_name,nrow=3) + 
  scale_fill_manual(values = color_vec_sbs, limits=force)  + 
  guides(fill = "none") + 
  labs(y="Proportional mutation source effect size",x="Signature")  + 
  theme(text = element_text(size = plot_text_size)) -> 
  fig1_attributable_ES



lusc_attributable_ES_long %>%
  group_by(bootstrap_sample,sig_label) %>%
  summarize(weight_prop = sum(weight_prop)) %>%
  group_by(sig_label) %>%
  summarize(avg_weight = mean(weight_prop)) -> 
  mean_attributable_effect


mean_attributable_effect %>%
  mutate(signature_process = case_when(
    sig_label == "SBS1" ~ "Deamination with age, clock-like (1)",
    sig_label == "SBS5" ~ "Unknown, clock-like (5)",
    sig_label %in% c("SBS2","SBS13") ~ "APOBEC (2,13)",
    sig_label %in% c("SBS3") ~ "Defective homologous recombination (3)",
    sig_label %in% c("SBS4","SBS29") ~ "Tobacco (4,29)",
    sig_label %in% c("SBS7a","SBS7b","SBS7c","SBS7d","SBS38") ~ "UV light (7a–d,38)",
    sig_label %in% c("SBS11","SBS31","SBS32","SBS35") ~ "Prior treatment (11,31,32,35)",
    sig_label %in% c("SBS22","SBS24","SBS42","SBS88") ~ "Mutagenic chemical exposure (22,24,42,88)",
    sig_label %in% c("SBS16") ~ "Alcohol-associated (16)",
    TRUE ~ "Non-actionable and unknown signatures")) %>%
  group_by(signature_process) %>%
  summarize(avg_weight = sum(avg_weight)) -> 
  mean_attributable_effect


mean_attributable_effect %>%
  arrange(avg_weight) %>%
  pull(signature_process) -> 
  sig_order_attribute

mean_attributable_effect$signature_process <- factor(mean_attributable_effect$signature_process, levels = sig_order_attribute)


ggplot(mean_attributable_effect, aes(fill=signature_process)) + 
  geom_bar(aes(x=1,y=avg_weight),stat="identity") + 
  scale_fill_manual(values=color_vec,limits=force) + 
  labs(y="Average attributable effect size\n among bootstrap samples",x=NULL,fill="Signature") + 
  theme_classic() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  guides(fill = "none") +
  theme(text = element_text(size = plot_text_size))-> 
  fig1_attributable_effectsize_stacked_nolegend


color_vec_fig1 <- color_vec[!names(color_vec) %in% c("Prior treatment (11,31,32,35)","Mutagenic chemical exposure (22,24,42,88)","Alcohol-associated (16)")]

ggplot(mean_attributable_effect, aes(fill=signature_process)) + 
  geom_bar(aes(x=1,y=avg_weight),stat="identity") + 
  scale_fill_manual(values=color_vec_fig1) + 
  labs(y="Average attributable effect size\n among bootstrap samples",x=NULL,fill="Signature") + 
  theme_classic() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.position = "bottom")  +
  theme(text = element_text(size = plot_text_size))-> 
  gg_for_legend

fig1_legend <- suppressWarnings(cowplot::get_legend(gg_for_legend))







