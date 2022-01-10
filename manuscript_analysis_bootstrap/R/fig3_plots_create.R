


# tumor_types <- names(fig3_data_gathered)
# 
# variant_attribution_plots <- vector(mode = "list",length = length(tumor_types))
# names(variant_attribution_plots) <- tumor_types

variant_attribution_plotter <- function(tumor_type, sig_of_focus=NULL){
  
  # tumor_type <- "LUAD" # for dev
  # sig_of_focus <- "Tobacco (4,29)"
  

  these_fig3_data <- fig3_data_gathered[[tumor_type]]
  
  
  these_fig3_medians <- these_fig3_data %>% 
    group_by(variant,signature) %>%
    summarize(median_weight = median(avg_weight)) %>%
    mutate(signature = str_remove(string = signature, pattern = "_nrsi"))
  
  
  these_fig3_medians_forplot <- these_fig3_medians %>%
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
    group_by(variant, signature_process) %>%
    summarize(median_weight = sum(median_weight)) %>% 
    ungroup() %>%
    mutate(value_to_plot = case_when(
      signature_process == sig_of_focus ~ median_weight*-1,
      TRUE ~ median_weight)) %>% 
    mutate(variant_name = variant_name_extracter(variant))
  
  
  these_fig3_medians_forplot %>%
    filter(!signature_process %in% sig_of_focus) %>%
    group_by(variant_name) %>%
    summarize(total = sum(median_weight)) %>%
    arrange(desc(total)) %>% 
    pull(variant_name) -> 
    order_of_variants
  
  these_fig3_medians_forplot$variant_name <- factor(these_fig3_medians_forplot$variant_name, 
                                                    levels = rev(order_of_variants))
  
  
  range_of_plot <- these_fig3_medians_forplot %>%
    mutate(collapsed_sig = case_when(
      signature_process == sig_of_focus ~ "below", 
      TRUE ~ "above")) %>% 
    group_by(variant_name, collapsed_sig) %>%
    summarize(total_bar = sum(median_weight))
  
  these_fig3_data %>% 
    mutate(signature = stringr::str_remove(string = signature, pattern = "_nrsi")) %>%
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
    filter(signature_process %in% sig_of_focus) %>% 
    mutate(variant_name = variant_name_extracter(variant)) %>%
    group_by(bootstrap_sample, variant_name, signature_process) %>%
    summarize(weight = sum(avg_weight)) %>% 
    group_by(variant_name) %>%
    summarize(qfirst = quantile(weight,0.025),
              qlast = quantile(weight,0.975)) -> 
    fig3_quantile_data
  
  
  these_fig3_data %>% 
    mutate(signature = stringr::str_remove(string = signature, pattern = "_nrsi")) %>%
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
    # filter(signature_process %in% sig_of_focus) %>% 
    mutate(variant_name = variant_name_extracter(variant)) %>%
    group_by(bootstrap_sample, variant_name, signature_process) %>%
    summarize(weight = sum(avg_weight)) %>% 
    group_by(variant_name, signature_process) %>%
    summarize(conf_lower = quantile(weight,0.025),
              median = median(weight),
              conf_upper = quantile(weight,0.975)) %>%
    mutate(tumor_type = tumor_type) -> 
    fig3_quantile_data_all
  
  
  
  plot_range <- max(max(range_of_plot$total_bar),fig3_quantile_data$qlast)
  
  y_axis_values <- pretty(c(-1*plot_range,plot_range),n = 4)
  
  if(length(y_axis_values) > 5){
    
    # only want text immediately around 0
    y_axis_values <- y_axis_values[y_axis_values %in% c(sort(abs(y_axis_values))[1:5],
                                                        sort(abs(y_axis_values))[1:5]*-1)]
    
  }
  
  these_fig3_medians_forplot$signature_process <- factor(these_fig3_medians_forplot$signature_process, levels = names(color_vec))
  
  this_fig3_attr_plot <- ggplot(these_fig3_medians_forplot, 
                                aes(y=variant_name)) + 
    geom_col(aes(x=value_to_plot,fill=signature_process),color="black") + 
    scale_fill_manual(values= color_vec, limits = force) + 
    guides(fill = "none") + 
    theme_bw() + 
    geom_vline(xintercept = 0) + 
    theme(axis.title = element_blank()) + 
    scale_x_continuous(breaks = y_axis_values,
                       labels = abs(y_axis_values),
                       limits = if(plot_range<y_axis_values[length(y_axis_values)]){
                         c(min(y_axis_values)-0.0,max(y_axis_values)+0.0)
                       }else{
                         c(-1*plot_range,plot_range)
                       }) + 
    theme(text = element_text(size = plot_text_size))
  
  
  this_fig3_attr_plot <- this_fig3_attr_plot + 
    geom_errorbar(data = fig3_quantile_data,aes(xmin=qfirst*-1,xmax=qlast*-1),width=0.5)
  
  
  
  return(list(this_fig3_attr_plot=this_fig3_attr_plot,
              fig3_quantile_data_all=fig3_quantile_data_all))
  
  
}



avg_signature_plotter <- function(tumor_type, sig_of_focus=NULL,plot_title=NULL){
  
  # ## Dev
  # tumor_type <- "LUAD"
  # sig_of_focus = "Tobacco (4,29)" 
  # plot_title <- "LUAD"
  # 
  
  these_trinucs <- fig3_trinuc_data_gathered[[tumor_type]]
  
  
  these_trinucs_forplot <- these_trinucs %>%
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
    group_by(signature_process) %>%
    summarize(avg_weight = sum(avg_weight)) %>%
    ungroup() %>%
    mutate(value_to_plot = case_when(
      signature_process == sig_of_focus ~ avg_weight*-1,
      TRUE ~ avg_weight))
  
  
  these_trinucs_forplot$signature_process <- factor(these_trinucs_forplot$signature_process, levels = names(color_vec))
  
  ggplot(these_trinucs_forplot) + 
    geom_col(aes(y=value_to_plot,x=1,fill=signature_process), color="black") + 
    coord_flip() + 
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits=c(-1,1),
                       # breaks = y_axis_values,
                       labels = c("1.0","0.5","0.0","0.5","1.0")) + 
    scale_fill_manual(values = color_vec,limits=force) + 
    guides(fill = "none") + 
    theme_bw() + 
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) + 
    labs(title = plot_title) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(text = element_text(size = plot_text_size))
  
  
}





# 
# variant_attribution_plotter(tumor_type = "LUAD",sig_of_focus = "Tobacco (4,29)")
# variant_attribution_plotter(tumor_type = "LUSC",sig_of_focus = "Tobacco (4,29)")
# variant_attribution_plotter(tumor_type = "SKCMP",sig_of_focus = "UV light (7a–d,38)")
# variant_attribution_plotter(tumor_type = "LIHC",sig_of_focus = "Mutagenic chemical exposure (22,24,42,88)")
# variant_attribution_plotter(tumor_type = "BLCA",sig_of_focus = "APOBEC (2,13)")
# variant_attribution_plotter(tumor_type = "CESC",sig_of_focus = "APOBEC (2,13)")
# variant_attribution_plotter(tumor_type = "HNSC_HPVneg",sig_of_focus = "APOBEC (2,13)")
# variant_attribution_plotter(tumor_type = "HNSC_HPVpos",sig_of_focus = "APOBEC (2,13)")
# 
# 
# 
# avg_signature_plotter(tumor_type = "HNSC_HPVpos",sig_of_focus = "APOBEC (2,13)",plot_title = "HNSC HPV negative")
# 




tumor_type <- "SKCMP" # for dev
sig_of_focus <- "UV light (7a–d,38)"

these_fig3_data <- fig3_data_gathered[[tumor_type]]


these_fig3_medians <- these_fig3_data %>% 
  group_by(variant,signature) %>%
  summarize(median_weight = median(avg_weight)) %>%
  mutate(signature = str_remove(string = signature, pattern = "_nrsi"))


these_fig3_medians_forplot <- these_fig3_medians %>%
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
  group_by(variant, signature_process) %>%
  summarize(median_weight = sum(median_weight)) %>% 
  ungroup() %>%
  mutate(value_to_plot = case_when(
    signature_process == sig_of_focus ~ median_weight*-1,
    TRUE ~ median_weight)) %>% 
  mutate(variant_name = variant_name_extracter(variant))


these_fig3_medians_forplot %>%
  filter(!signature_process %in% sig_of_focus) %>%
  group_by(variant_name) %>%
  summarize(total = sum(median_weight)) %>%
  arrange(desc(total)) %>% 
  pull(variant_name) -> 
  order_of_variants

these_fig3_medians_forplot$variant_name <- factor(these_fig3_medians_forplot$variant_name, 
                                                  levels = rev(order_of_variants))
these_fig3_medians_forplot

these_fig3_medians_forplot %>% 
  filter(variant_name == "KIT K642E") %>%
  mutate(med_weight_prop = median_weight / sum(median_weight)) %>% 
  select(-value_to_plot) %>% arrange(desc(med_weight_prop)) %>%
  select(signature_process, med_weight_prop)




