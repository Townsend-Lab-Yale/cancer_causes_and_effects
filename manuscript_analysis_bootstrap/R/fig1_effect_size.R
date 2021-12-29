# combine effect size data for that specific tumor

lusc_ces_rbind <- rbindlist(lusc_ces,idcol = "bootstrap_sample")


lusc_ces_rbind %>%
  select(bootstrap_sample,variant_id,selection_intensity) %>%
  mutate(variant_name = variant_name_extracter(variant_id)) -> 
  effect_size_long 


effect_size_long %>%
  group_by(variant_name) %>%
  summarize(qfirst = quantile(selection_intensity,0.025),
                        median = median(selection_intensity),
                        mean = mean(selection_intensity),
                        qlast = quantile(selection_intensity, 0.975)) -> 
  effect_size_summary

effect_size_summary$variant_name <- factor(effect_size_summary$variant_name,levels = c("TP53 R282W",
                                                                                       "KLF5 E419Q",
                                                                                       "OR2T34 L163L"))


ggplot(effect_size_summary, aes(x=variant_name, y= median,fill = variant_name)) + 
  geom_col() + 
  geom_errorbar(aes(ymin=qfirst,ymax=qlast),width=0.3) + 
  geom_jitter(data = effect_size_long, aes(y=selection_intensity),shape = 21,alpha=0.2,height = 0) + 
  scale_fill_manual(values = color_vec,limits = force) + 
  theme_classic() + 
  guides(fill="none") + 
  labs(x="Variant",y="Effect size") +
  theme(text = element_text(size = plot_text_size),
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) -> 
  fig1_effect_size






