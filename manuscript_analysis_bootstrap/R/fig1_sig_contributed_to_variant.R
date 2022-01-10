

# script to visualize the probability each source contributed to each tumor





variant_name_extracter <- function(variant_names_vec){
  variant_names <- strsplit(x = variant_names_vec, split = "_")
  variant_names <- unlist(lapply(variant_names, function(x) paste(c(x[1],x[2]),collapse = " ")))
  return(variant_names)
}

lusc_attribution_long <- lusc_attribution_long %>%
  mutate(variant_name = variant_name_extracter(variant))

lusc_attribution_long$variant_name <- factor(lusc_attribution_long$variant_name, levels = c("TP53 R282W",
                                                                                            "KLF5 E419Q",
                                                                                            "OR2T34 L163L"))

lusc_attribution_long <- lusc_attribution_long %>%
  mutate(sig_label = stringr::str_remove(string = signature, pattern = "_flux"))

lusc_attribution_long$sig_label <- factor(lusc_attribution_long$sig_label, levels = sig_order)


lusc_attribution_long %>%
  group_by(sig_label,variant_name) %>%
  summarize(qfirst = quantile(weight,0.025),
            median = median(weight),
            mean = mean(weight),
            qlast = quantile(weight, 0.975)) -> 
  contribution_summary

contribution_summary %>%
  filter(qlast>0) %>%
  pull(sig_label) %>%
  as.character() -> 
  sigs_to_plot

lusc_attribution_long <- lusc_attribution_long %>%
  filter(sig_label %in% sigs_to_plot)


# lusc_attribution_long %>% 
#   ggplot(aes(x=sig_label,y=weight,fill = sig_label)) + 
#   geom_jitter(shape = 21, alpha=0.2) +
#   scale_fill_manual(values  = color_vec_sbs,limits = force) + 
#   facet_wrap(~variant_name) + 
#   theme_classic() + 
#   theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) + 
#   guides(fill="none")

contribution_summary <- contribution_summary %>%
  filter(sig_label %in%  sigs_to_plot)

x_labels <- as.character(sort(unique(contribution_summary$sig_label))) %>%
  stringr::str_remove(pattern = "SBS")

borders_colors <- color_vec[levels(contribution_summary$variant_name)]


ggplot(contribution_summary, aes(x=sig_label,fill=sig_label)) + 
  geom_col(aes(y=median),color="black") + 
  geom_jitter(data = lusc_attribution_long, aes(y=weight), shape = 21, alpha = 0.2, width = 0.2,height = 0) + 
  geom_errorbar(aes(ymin = qfirst, ymax = qlast), width = 0.3) + 
  theme_classic() + 
  scale_x_discrete(labels = x_labels) + 
  facet_wrap(~variant_name,nrow=3) + 
  scale_fill_manual(values = color_vec_sbs, limits=force)  + 
  guides(fill = "none") + 
  theme(text = element_text(size = plot_text_size)) +
  labs(y="Probability each source\n contributed to each variant",x="Signature")  -> 
  fig1_contributions






# 
# ggplot(weight_summary, aes(x=signature,fill=signature)) + 
#   geom_col(aes(y=median),color="black") + 
#   geom_jitter(data = lusc_tumor_trinuc_long, aes(y=weight),shape=21,alpha=0.2,width = 0.2) +
#   geom_errorbar(aes(ymin = qfirst,ymax=qlast), width = 0.3) + 
#   theme_classic() + 
#   scale_x_discrete(labels = x_labels) + 
#   # limits = force gets rid of extra colors
#   scale_fill_manual(values = color_vec_sbs, limits=force) + 
#   guides(fill = "none") + 
#   labs(y="Signature weight",x="Signature") -> 
#   fig1_signatures_in_tumor
#                                  
                                 
                                 
                                 
                                 



