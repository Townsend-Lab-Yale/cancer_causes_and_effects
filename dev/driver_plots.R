# 
# 
# 
# # colors
# # color_vec_fullname
# 

library(tidyverse)
library(cancereffectsizeR)

# load in data ----
load("input_data/from_cluster/MP_2_2_0-all_scaled_selection_main.RData")
load("input_data/from_cluster/MP_2_2_0-combined_selection_results.RData")

# the following file is from this manuscript:
# Bailey, Matthew H., Collin Tokheim, Eduard Porta-Pardo, Sohini Sengupta, Denis Bertrand, Amila Weerasinghe, Antonio Colaprico, et al. 2018. “Comprehensive Characterization of Cancer Driver Genes and Mutations.” Cell 174 (4): 1034–35.

# The file is located here: https://www.cell.com/cms/10.1016/j.cell.2018.02.060/attachment/b0abfafa-a1ff-4f7d-8842-b49ba8d32e08/mmc1.xlsx

# I will not include the file in our repo as I am unsure if we have permission to
# share/host. But, it is open access!
bailey_driver_list <- readxl::read_excel(path = "dev/Bailey_etal_Cell_1-s2.0-S009286741830237X-mmc1.xlsx",sheet = "Table S1",skip = 3)


load("input_data/from_cluster/MP_2_2_0-weights_combined.RData")
load("input_data/from_cluster/MP_2_2_0-variant_prevalence_main.RData")



# setting up data ----

# scaled selection
all_scaled_selection_main_m <- pivot_longer(data = all_scaled_selection_main,
                                            cols = starts_with("SBS"))

all_scaled_selection_main_m <- all_scaled_selection_main_m %>%
  filter(value>0) %>%
  separate(col = name, into = c("Signature","Type"),sep = "_")

all_scaled_selection_main_m <- all_scaled_selection_main_m %>%
  mutate(variant = stringr::str_replace(string = variant, pattern = "_ENSP[:digit:]+",replacement = ""))


# weights


# sum(weights_combined$Weight_prop[weights_combined$Tumor=="TCGA-2F-A9KO"])

avg_weights <- weights_combined %>%
  group_by(tumor_type,Signature) %>%
  summarize(avg_weight = mean(Weight_prop))

signatures_names_matrix <- cancereffectsizeR::get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.1")

avg_weights_filtered <- avg_weights %>%
  filter(avg_weight >0)

avg_weights_filtered$Signature <- factor(avg_weights_filtered$Signature,levels = signatures_names_matrix$meta$Signature[
  signatures_names_matrix$meta$Signature %in% avg_weights_filtered$Signature])

avg_weights_filtered$Sig_num <- gsub(pattern = "SBS",replacement = "",as.character(avg_weights_filtered$Signature))



# variant specific nrsi----
scaled_selection <- all_scaled_selection_main_m %>%
  filter(Type == "nrsi")


scaled_selection <- scaled_selection %>%
  mutate(variant_original = variant) %>%
  separate(variant,into = c("Gene","AA_change"),sep = "_")

scaled_selection <- scaled_selection %>%
  filter(!is.na(value))


scaled_selection <- scaled_selection %>%
  group_by(Unique_Patient_Identifier) %>%
  mutate(nrsi_per_tumor = value/sum(value))











signatures_names_matrix <- cancereffectsizeR::get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.1")
signatures_names_matrix <- as.data.frame(signatures_names_matrix$meta[,c("Signature","Etiology")])

signatures_names_matrix$Etiology_sig <- paste0(signatures_names_matrix$Etiology, " (",gsub(pattern = "SBS",replacement = "",x = signatures_names_matrix$Signature),")")

rownames(signatures_names_matrix) <- signatures_names_matrix$Signature


signatures_names_matrix$sig_num <- gsub(pattern = "SBS",replacement = "",x = signatures_names_matrix$Signature)





library(tidyverse)
library(cancereffectsizeR) # use v2.0
library(cowplot)
library(RColorBrewer)

variant_name_extracter <- function(variant_names_vec){
  variant_names <- strsplit(x = variant_names_vec, split = "_")
  variant_names <- unlist(lapply(variant_names, function(x) paste(c(x[1],x[2]),collapse = " ")))
  return(variant_names)
}

color_vec_fullname <- setNames(object = c("gray60","gray80",RColorBrewer::brewer.pal(n = 7,name = "Set3"),"black")
                               , nm = c("Deamination with age, clock-like (1)",
                                        "Unknown, clock-like (5)",
                                        "APOBEC (2,13)",
                                        "Defective homologous recombination (3)",
                                        "Tobacco (4,29)",
                                        "UV light (7a–d,38)",
                                        "Prior treatment (11,31,32,35)",
                                        "Mutagenic chemical exposure (22,24,42,88)",
                                        "Alcohol-associated (16)",
                                        "Non-actionable or unknown signatures"
                               ))

# switching two so UV light is the bright yellow
color_vec_fullname["UV light (7a–d,38)"] <- "#FFFFB3"
color_vec_fullname["Defective homologous recombination (3)"] <- "#FB8072"






# # function to plot driver focused nrsi---- 
# driver_focused_nrsi_weight_barplot <- function(tumor_type_ourdata,
#                                                Bailey_tumor_type,
#                                                signatures_below_line,
#                                                signatures_below_line_grouped,
#                                                drivers_to_plot=10,
#                                                ordered_by_total_weight=F,
#                                                text_font_size=10,
#                                                legend_columns = 2,
#                                                below_line_black=F,
#                                                choose_variants_by_total_volume=F,
#                                                plot_title=NULL){
#   
#   scaled_selection %>%
#     filter(tumor_type == tumor_type_ourdata) %>%
#     filter(Gene %in% bailey_driver_list$Gene[bailey_driver_list$Cancer==Bailey_tumor_type]) %>%
#     # filter(variant_original == "BRAF_V600E") %>%
#     mutate(Signature = forcats::as_factor(Signature)) %>%
#     # mutate(sig_collapsed = forcats::fct_collapse(Signature,
#     #                                              SBS7_all = c("SBS7a","SBS7b","SBS7c","SBS7d"))) %>%
#     # mutate(sig_collapsed = forcats::fct_collapse(sig_collapsed,
#     #                                              `SBS2,13` = c("SBS2","SBS13"))) %>%
#     # tidyr::complete(sig_collapsed, fill = list(nrsi_per_tumor = 0)) %>%
#     group_by(Unique_Patient_Identifier,variant_original,Signature) %>%
#     summarize(nrsi_per_tumor = sum(nrsi_per_tumor)) %>%
#     ungroup() %>%
#     group_by(variant_original,Signature) %>%
#     summarize(summed_sig_nrsi = sum(nrsi_per_tumor)) %>%
#     ungroup() -> driver_signature_nrsi_weight
#   
#   
#   if(choose_variants_by_total_volume){
#     
#     drivers_to_pick <- driver_signature_nrsi_weight %>%
#       group_by(variant_original) %>%
#       summarize(total_volume = sum(summed_sig_nrsi)) %>%
#       arrange(desc(total_volume)) %>%
#       pull(variant_original) %>%
#       .[1:drivers_to_plot]
#     
#   }else{
#     drivers_to_pick <- driver_signature_nrsi_weight$variant_original
#   }
#   
#   if(ordered_by_total_weight){
#     order_of_drivers <- driver_signature_nrsi_weight %>%
#       filter(variant_original %in% drivers_to_pick) %>%
#       group_by(variant_original) %>%
#       summarize(summed_sigs = sum(summed_sig_nrsi)) %>%
#       arrange(desc(summed_sigs)) %>%
#       pull(variant_original) %>%
#       .[!is.na(.)] %>%
#       .[1:drivers_to_plot]
#   }else{
#     order_of_drivers <- driver_signature_nrsi_weight %>%
#       filter(variant_original %in% drivers_to_pick) %>%
#       filter(!Signature %in% signatures_below_line) %>%
#       group_by(variant_original) %>%
#       summarize(summed_sigs = sum(summed_sig_nrsi)) %>%
#       arrange(desc(summed_sigs)) %>%
#       pull(variant_original) %>%
#       .[!is.na(.)] %>%
#       .[1:drivers_to_plot]
#   }
#   
#   
#   
#   driver_signature_nrsi_weight %>%
#     mutate(signature_full = signatures_names_matrix[as.character(Signature), "Etiology_sig"]) %>%
#     mutate(signature_full =
#              forcats::fct_collapse(
#                signature_full,
#                `UV light (7a–d,38)` = c("UV light (7a)",
#                                         "UV light (7b)",
#                                         "UV light (7c)",
#                                         "UV light (7d)",
#                                         "Potentially indirect damage from UV light (38)"),
#                `APOBEC (2,13)` = c("APOBEC (2)","APOBEC (13)"),
#                `Tobacco (4,29)` = c("Tobacco smoking (4)", "Tobacco chewing (29)"),
#                `Prior treatment (11,31,32,35)` = c("Temozolomide treatment (11)",
#                                                    "Platinum drug chemotherapy (31)",
#                                                    "Azathioprine treatment (used for immunosuppression) (32)",
#                                                    "Prior chemotherapy treatment (35)"),
#                `Mutagenic chemical exposure (22,24,42,88)` = c("Aristolochic acid exposure (22)",
#                                                                "Aflatoxin exposure (24)",
#                                                                "Occupational exposure to haloalkanes (42)",
#                                                                "Colibactin exposure (COSMIC 3.1) (88)")
#                
#              )
#     ) %>%
#     mutate(signature_full =
#              forcats::fct_recode(signature_full,
#                                  `Deamination with age, clock-like (1)` = "Deamination with age (1)",
#                                  `Alcohol-associated (16)` = "Unknown (16)")) ->
#     driver_signature_nrsi_weight
#   
#   
#   driver_signature_nrsi_weight %>%
#     left_join(
#       filter(
#         variant_prevalence_main,
#         tumor_type == tumor_type_ourdata),
#       by = c("variant_original" = "variant_name")) ->
#     driver_signature_nrsi_weight
#   
#   
#   
#   
#   
#   # mutate(signature_full = case_when(
#   #   as.character(signature_full) %in% names(color_vec_fullname) ~ signature_full,
#   #   TRUE ~ "Non-actionable or unknown signatures")) %>%
#   driver_signature_nrsi_weight %>%
#     group_by(variant_original,signature_full) %>%
#     summarize(avg_sig_nrsi = sum(nrsi_per_tumor)/samples_covering) %>%
#     ungroup() ->
#     driver_signature_nrsi_weight_full_sig
#   
#   
#   unique(driver_signature_nrsi_weight_full_sig$signature_full)
#   
#   driver_signature_nrsi_weight_forplot <- driver_signature_nrsi_weight_full_sig %>%
#     filter(variant_original %in% order_of_drivers) %>%
#     mutate(variant_original = factor(variant_original,levels=order_of_drivers)) %>%
#     mutate(avg_sig_nrsi = case_when(signature_full %in% signatures_below_line_grouped ~ -avg_sig_nrsi,
#                                     TRUE ~ avg_sig_nrsi))
#   
#   
#   driver_signature_nrsi_weight_forplot$signature_full <- forcats::fct_relevel(
#     driver_signature_nrsi_weight_forplot$signature_full,
#     names(color_vec_fullname))
#   
#   # y_axis_values <- pretty(driver_signature_nrsi_weight_forplot$summed_sig_nrsi)
#   
#   
#   driver_signature_nrsi_weight_forplot %>%
#     mutate(sig_char = as.character(signature_full)) %>%
#     mutate(sig_char = case_when(
#       sig_char %in% names(color_vec_fullname) ~ sig_char,
#       TRUE ~ "Non-actionable or unknown signatures")) %>%
#     mutate(sig_char = forcats::fct_relevel(sig_char,
#                                            names(color_vec_fullname))) ->
#     driver_signature_nrsi_weight_forplot
#   
#   
#   nrsi_summed <- driver_signature_nrsi_weight_forplot %>%
#     mutate(collapsed_sig = case_when(
#       sig_char == signatures_below_line_grouped ~ "below" ,
#       TRUE ~ "above")) %>%
#     group_by(variant_original,collapsed_sig) %>%
#     summarize(main_bar = sum(abs(avg_sig_nrsi)))
#   
#   
#   plot_range <- max(nrsi_summed$main_bar)
#   
#   y_axis_values <- pretty(c(-1*plot_range,plot_range))
#   
#   
#   driver_signature_nrsi_weight_forplot$variant_original <-
#     driver_signature_nrsi_weight_forplot$variant_original %>%
#     fct_relabel(~ gsub(pattern = "_",replacement = " ", .x))
#   
#   
#   driver_plot_out <- ggplot(driver_signature_nrsi_weight_forplot) +
#     geom_bar(aes(x=variant_original,y = avg_sig_nrsi,fill=sig_char),stat="identity") +
#     scale_fill_manual(values = color_vec_fullname, na.value="black") +
#     theme_bw() +
#     geom_hline(yintercept = 0) +
#     scale_y_continuous(breaks = y_axis_values,
#                        labels = abs(y_axis_values),
#                        limits = c(-1*plot_range,plot_range))+
#     theme(axis.text.x = element_text(angle=35,vjust=1,hjust=1)) +
#     labs(title=
#            paste(tumor_type_ourdata),
#          y="Summed per-tumor relative variant effect size\ncontributed from each signature weight",
#          x="Variant",
#          fill = "Signature") +
#     theme(text = element_text(size=text_font_size)) +
#     guides(fill=guide_legend(ncol=legend_columns))
#   
#   
#   # flipped
#   
#   driver_plot_out <- ggplot(driver_signature_nrsi_weight_forplot) +
#     geom_bar(aes(y=fct_rev(variant_original),x = avg_sig_nrsi,fill=sig_char),stat="identity") +
#     scale_fill_manual(values = color_vec_fullname, na.value="black") +
#     theme_bw() +
#     geom_vline(xintercept = 0) +
#     scale_x_continuous(breaks = y_axis_values,
#                        labels = abs(y_axis_values),
#                        limits = c(-1*plot_range,plot_range))+
#     theme(axis.text.y = element_text(angle=0,vjust=0.5,hjust=1)) +
#     # scale_y_discrete(position = "right") +
#     labs(
#       x="Summed per-tumor relative variant effect size\ncontributed from each signature weight",
#       fill = "Signature") +
#     theme(text = element_text(size=text_font_size)) +
#     guides(fill=guide_legend(ncol=legend_columns))
#   
#   
#   # driver_plot_out + theme(legend.position = "none")
#   
#   # avg weight plot ----
#   
#   
#   avg_weights %>%
#     ungroup() %>%
#     filter(tumor_type == tumor_type_ourdata) %>%
#     mutate(sig_full = signatures_names_matrix[as.character(Signature),"Etiology_sig"]) %>%
#     mutate(sig_full_fct =
#              forcats::fct_collapse(sig_full,
#                                    `UV light (7a–d,38)` = c("UV light (7a)",
#                                                             "UV light (7b)",
#                                                             "UV light (7c)",
#                                                             "UV light (7d)",
#                                                             "Potentially indirect damage from UV light (38)"),
#                                    `APOBEC (2,13)` = c("APOBEC (2)","APOBEC (13)"),
#                                    `Tobacco (4,29)` = c("Tobacco smoking (4)", "Tobacco chewing (29)"),
#                                    `Prior treatment (11,31,32,35)` = c("Temozolomide treatment (11)",
#                                                                        "Platinum drug chemotherapy (31)",
#                                                                        "Azathioprine treatment (used for immunosuppression) (32)",
#                                                                        "Prior chemotherapy treatment (35)"),
#                                    `Mutagenic chemical exposure (22,24,42,88)` = c("Aristolochic acid exposure (22)",
#                                                                                    "Aflatoxin exposure (24)",
#                                                                                    "Occupational exposure to haloalkanes (42)",
#                                                                                    "Colibactin exposure (COSMIC 3.1) (88)")
#                                    
#              )
#     ) %>%
#     mutate(sig_full_fct =
#              forcats::fct_recode(sig_full_fct,
#                                  `Deamination with age, clock-like (1)` = "Deamination with age (1)",
#                                  `Alcohol-associated (16)` = "Unknown (16)")) %>%
#     group_by(tumor_type,sig_full_fct) %>%
#     summarize(avg_weight = sum(avg_weight)) ->
#     avg_weights_sig_full
#   
#   
#   avg_weights_sig_full <- avg_weights_sig_full %>%
#     mutate(avg_weight = case_when(sig_full_fct %in% signatures_below_line_grouped ~ -avg_weight,
#                                   TRUE ~ avg_weight))
#   
#   
#   avg_weights_sig_full %>%
#     mutate(sig_char = as.character(sig_full_fct)) %>%
#     mutate(sig_char = case_when(
#       sig_char %in% names(color_vec_fullname) ~ sig_char,
#       TRUE ~ "Non-actionable or unknown signatures")) %>%
#     mutate(sig_char = forcats::fct_relevel(sig_char,
#                                            names(color_vec_fullname))) ->
#     avg_weights_sig_full_forplot
#   
#   
#   avg_weight_plot_out <- ggplot(avg_weights_sig_full_forplot) +
#     geom_bar(aes(x=1,y = avg_weight,fill=sig_char),stat="identity") +
#     scale_fill_manual(values = color_vec_fullname, na.value="black") +
#     theme_bw() +
#     geom_hline(yintercept = 0) +
#     scale_y_continuous(limits=c(-1,1)
#                        # breaks = y_axis_values,
#                        # labels = abs(y_axis_values)
#     )+
#     theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
#     labs(title= plot_title,
#          x="Average signature weight",
#          # x="Variant",
#          fill = "Signature") +
#     theme(text = element_text(size=text_font_size),
#           axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
#     guides(fill=guide_legend(ncol=legend_columns)) +
#     theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
#   
#   
#   
#   
#   avg_weight_plot_out <- ggplot(avg_weights_sig_full_forplot) +
#     geom_bar(aes(x=1,y = avg_weight,fill=sig_char),stat="identity") +
#     coord_flip() +
#     scale_fill_manual(values = color_vec_fullname, na.value="black") +
#     theme_bw() +
#     geom_hline(yintercept = 0) +
#     scale_y_continuous(limits=c(-1,1)
#                        # breaks = y_axis_values,
#                        # labels = abs(y_axis_values)
#     )+
#     theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=1)) +
#     labs(title = plot_title,
#          y="Average signature weight",
#          # x="Variant",
#          fill = "Signature") +
#     theme(text = element_text(size=text_font_size),
#           axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
#     guides(fill=guide_legend(ncol=legend_columns)) +
#     theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
#   
#   
#   
#   return(list(effect_plot = driver_plot_out + theme(legend.position = "none"),
#               avg_weight_plot = avg_weight_plot_out+ theme(legend.position = "none")))
#   
# }
# 
# 
# 
# 
# 
# 
# 

# function to plot driver focused nrsi---- 
driver_focused_nrsi_weight_barplot_avg <- function(tumor_type_ourdata,
                                                   Bailey_tumor_type,
                                                   signatures_below_line,
                                                   signatures_below_line_grouped,
                                                   drivers_to_plot=10,
                                                   ordered_by_total_weight=F,
                                                   text_font_size=10,
                                                   legend_columns = 2,
                                                   below_line_black=F,
                                                   choose_variants_by_total_volume=F,
                                                   plot_title=NULL){
  
  scaled_selection %>%
    filter(tumor_type == tumor_type_ourdata) %>%
    filter(Gene %in% bailey_driver_list$Gene[bailey_driver_list$Cancer==Bailey_tumor_type]) %>%
    # filter(variant_original == "BRAF_V600E") %>%
    mutate(Signature = forcats::as_factor(Signature)) -> driver_signature_nrsi_weight
  
  
  # %>%
  #   # mutate(sig_collapsed = forcats::fct_collapse(Signature, 
  #   #                                              SBS7_all = c("SBS7a","SBS7b","SBS7c","SBS7d"))) %>% 
  #   # mutate(sig_collapsed = forcats::fct_collapse(sig_collapsed, 
  #   #                                              `SBS2,13` = c("SBS2","SBS13"))) %>% 
  #   # tidyr::complete(sig_collapsed, fill = list(nrsi_per_tumor = 0)) %>%
  #   group_by(Unique_Patient_Identifier,variant_original,Signature) %>%
  #   summarize(nrsi_per_tumor = sum(nrsi_per_tumor)) %>%
  #   ungroup() %>%
  #   group_by(variant_original,Signature) %>%
  #   summarize(summed_sig_nrsi = sum(nrsi_per_tumor)) %>%
  #   ungroup() -> driver_signature_nrsi_weight
  
  
  if(choose_variants_by_total_volume){
    
    drivers_to_pick <- driver_signature_nrsi_weight %>%
      group_by(variant_original) %>%
      summarize(total_volume = sum(nrsi_per_tumor)) %>%
      arrange(desc(total_volume)) %>%
      pull(variant_original) %>%
      .[1:drivers_to_plot]
    
  }else{
    drivers_to_pick <- driver_signature_nrsi_weight$variant_original
  }
  
  if(ordered_by_total_weight){
    order_of_drivers <- driver_signature_nrsi_weight %>%
      filter(variant_original %in% drivers_to_pick) %>%
      group_by(variant_original) %>%
      summarize(summed_sigs = sum(nrsi_per_tumor)) %>%
      arrange(desc(summed_sigs)) %>%
      pull(variant_original) %>%
      .[!is.na(.)] %>%
      .[1:drivers_to_plot]
  }else{
    order_of_drivers <- driver_signature_nrsi_weight %>%
      filter(variant_original %in% drivers_to_pick) %>%
      filter(!Signature %in% signatures_below_line) %>%
      group_by(variant_original) %>%
      summarize(summed_sigs = sum(nrsi_per_tumor)) %>%
      arrange(desc(summed_sigs)) %>%
      pull(variant_original) %>%
      .[!is.na(.)] %>%
      .[1:drivers_to_plot]
    
    if(any(is.na(order_of_drivers))){
     
      re_include <- setdiff(drivers_to_pick,order_of_drivers)
      
      
      if(length(re_include)>0){
        order_of_drivers[is.na(order_of_drivers)] <- re_include
        }
      
    }
  }
  
  
  driver_signature_nrsi_weight %>%
    mutate(signature_full = signatures_names_matrix[as.character(Signature), "Etiology_sig"]) %>%
    mutate(signature_full =
             suppressWarnings(forcats::fct_collapse(
               signature_full,
               `UV light (7a–d,38)` = c("UV light (7a)",
                                        "UV light (7b)",
                                        "UV light (7c)",
                                        "UV light (7d)",
                                        "Potentially indirect damage from UV light (38)"),
               `APOBEC (2,13)` = c("APOBEC (2)","APOBEC (13)"),
               `Tobacco (4,29)` = c("Tobacco smoking (4)", "Tobacco chewing (29)"),
               `Prior treatment (11,31,32,35)` = c("Temozolomide treatment (11)",
                                                   "Platinum drug chemotherapy (31)",
                                                   "Azathioprine treatment (used for immunosuppression) (32)",
                                                   "Prior chemotherapy treatment (35)"),
               `Mutagenic chemical exposure (22,24,42,88)` = c("Aristolochic acid exposure (22)",
                                                               "Aflatoxin exposure (24)",
                                                               "Occupational exposure to haloalkanes (42)",
                                                               "Colibactin exposure (COSMIC 3.1) (88)")
               
             ))
    ) %>%
    mutate(signature_full =
             suppressWarnings(forcats::fct_recode(signature_full,
                                 `Deamination with age, clock-like (1)` = "Deamination with age (1)",
                                 `Alcohol-associated (16)` = "Unknown (16)"))) ->
    driver_signature_nrsi_weight
  
  
  driver_signature_nrsi_weight %>%
    left_join(
      filter(
        variant_prevalence_main,
        tumor_type == tumor_type_ourdata),
      by = c("variant_original" = "variant_name")) ->
    driver_signature_nrsi_weight
  
  
  
  
  
  # mutate(signature_full = case_when(
  #   as.character(signature_full) %in% names(color_vec_fullname) ~ signature_full,
  #   TRUE ~ "Non-actionable or unknown signatures")) %>%
  driver_signature_nrsi_weight %>%
    group_by(variant_original,signature_full,samples_covering) %>%
    summarize(summed_sig_nrsi = sum(nrsi_per_tumor)) %>%
    ungroup() %>%
    mutate(avg_sig_nrsi = summed_sig_nrsi/samples_covering)->
    driver_signature_nrsi_weight_full_sig
  
  
  unique(driver_signature_nrsi_weight_full_sig$signature_full)
  
  driver_signature_nrsi_weight_forplot <- driver_signature_nrsi_weight_full_sig %>%
    filter(variant_original %in% order_of_drivers) %>%
    mutate(variant_original = factor(variant_original,levels=order_of_drivers)) %>%
    mutate(avg_sig_nrsi = case_when(signature_full %in% signatures_below_line_grouped ~ -avg_sig_nrsi,
                                    TRUE ~ avg_sig_nrsi))
  
  
  driver_signature_nrsi_weight_forplot$signature_full <- forcats::fct_relevel(
    driver_signature_nrsi_weight_forplot$signature_full,
    names(color_vec_fullname))
  
  # y_axis_values <- pretty(driver_signature_nrsi_weight_forplot$summed_sig_nrsi)
  
  
  driver_signature_nrsi_weight_forplot %>%
    mutate(sig_char = as.character(signature_full)) %>%
    mutate(sig_char = case_when(
      sig_char %in% names(color_vec_fullname) ~ sig_char,
      TRUE ~ "Non-actionable or unknown signatures")) %>%
    mutate(sig_char = forcats::fct_relevel(sig_char,
                                           names(color_vec_fullname))) ->
    driver_signature_nrsi_weight_forplot
  
  
  nrsi_summed <- driver_signature_nrsi_weight_forplot %>%
    mutate(collapsed_sig = case_when(
      sig_char == signatures_below_line_grouped ~ "below" ,
      TRUE ~ "above")) %>%
    group_by(variant_original,collapsed_sig) %>%
    summarize(main_bar = sum(abs(avg_sig_nrsi)))
  
  
  plot_range <- max(nrsi_summed$main_bar)
  
  y_axis_values <- pretty(c(-1*plot_range,plot_range),n = 3)
  
  if(length(y_axis_values) > 3){
    
    # only want text immediately around 0
    y_axis_values <- y_axis_values[y_axis_values %in% c(sort(abs(y_axis_values))[1:3],
                                                        sort(abs(y_axis_values))[1:3]*-1)]
    
  }
  
  driver_signature_nrsi_weight_forplot$variant_original <-
    driver_signature_nrsi_weight_forplot$variant_original %>%
    fct_relabel(~ gsub(pattern = "_",replacement = " ", .x))
  
  
  driver_plot_out <- ggplot(driver_signature_nrsi_weight_forplot) +
    geom_bar(aes(x=variant_original,y = avg_sig_nrsi,fill=sig_char),stat="identity") +
    scale_fill_manual(values = color_vec_fullname, na.value="black") +
    theme_bw() +
    geom_hline(yintercept = 0) +
    scale_y_continuous(breaks = y_axis_values,
                       labels = abs(y_axis_values),
                       limits = c(-1*plot_range,plot_range))+
    theme(axis.text.x = element_text(angle=35,vjust=1,hjust=1)) +
    labs(title=
           paste(tumor_type_ourdata),
         y="Summed per-tumor relative variant effect size\ncontributed from each signature weight",
         x="Variant",
         fill = "Signature") +
    theme(text = element_text(size=text_font_size)) +
    guides(fill=guide_legend(ncol=legend_columns))
  
  
  # flipped
  
  driver_plot_out <- ggplot(driver_signature_nrsi_weight_forplot) +
    geom_bar(aes(y=fct_rev(variant_original),
                 x = avg_sig_nrsi,fill=sig_char),
             stat="identity", color="gray20",size=0.1) +
    scale_fill_manual(values = color_vec_fullname, na.value="black") +
    theme_bw() +
    geom_vline(xintercept = 0) +
    scale_x_continuous(breaks = y_axis_values,
                       labels = abs(y_axis_values),
                       limits = c(-1*plot_range,plot_range))+
    theme(axis.text.y = element_text(angle=0,vjust=0.5,hjust=1)) +
    # scale_y_discrete(position = "right") +
    labs(
      x="Average per-tumor relative variant effect size\ncontributed from each signature weight",
      fill = "Signature") +
    theme(text = element_text(size=text_font_size)) +
    guides(fill=guide_legend(ncol=legend_columns))
  
  
  # driver_plot_out + theme(legend.position = "none")
  
  # avg weight plot ----
  
  
  avg_weights %>%
    ungroup() %>%
    filter(tumor_type == tumor_type_ourdata) %>%
    mutate(sig_full = signatures_names_matrix[as.character(Signature),"Etiology_sig"]) %>%
    mutate(sig_full_fct =
             suppressWarnings( forcats::fct_collapse(sig_full,
                                   `UV light (7a–d,38)` = c("UV light (7a)",
                                                            "UV light (7b)",
                                                            "UV light (7c)",
                                                            "UV light (7d)",
                                                            "Potentially indirect damage from UV light (38)"),
                                   `APOBEC (2,13)` = c("APOBEC (2)","APOBEC (13)"),
                                   `Tobacco (4,29)` = c("Tobacco smoking (4)", "Tobacco chewing (29)"),
                                   `Prior treatment (11,31,32,35)` = c("Temozolomide treatment (11)",
                                                                       "Platinum drug chemotherapy (31)",
                                                                       "Azathioprine treatment (used for immunosuppression) (32)",
                                                                       "Prior chemotherapy treatment (35)"),
                                   `Mutagenic chemical exposure (22,24,42,88)` = c("Aristolochic acid exposure (22)",
                                                                                   "Aflatoxin exposure (24)",
                                                                                   "Occupational exposure to haloalkanes (42)",
                                                                                   "Colibactin exposure (COSMIC 3.1) (88)")
                                   
             ) )
    ) %>%
    mutate(sig_full_fct =
             suppressWarnings(forcats::fct_recode(sig_full_fct,
                                 `Deamination with age, clock-like (1)` = "Deamination with age (1)",
                                 `Alcohol-associated (16)` = "Unknown (16)"))) %>%
    group_by(tumor_type,sig_full_fct) %>%
    summarize(avg_weight = sum(avg_weight)) ->
    avg_weights_sig_full
  
  
  avg_weights_sig_full <- avg_weights_sig_full %>%
    mutate(avg_weight = case_when(sig_full_fct %in% signatures_below_line_grouped ~ -avg_weight,
                                  TRUE ~ avg_weight))
  
  
  avg_weights_sig_full %>%
    mutate(sig_char = as.character(sig_full_fct)) %>%
    mutate(sig_char = case_when(
      sig_char %in% names(color_vec_fullname) ~ sig_char,
      TRUE ~ "Non-actionable or unknown signatures")) %>%
    mutate(sig_char = forcats::fct_relevel(sig_char,
                                           names(color_vec_fullname))) ->
    avg_weights_sig_full_forplot
  
  
  # avg_weight_plot_out <- ggplot(avg_weights_sig_full_forplot) +
  #   geom_bar(aes(x=1,y = avg_weight,fill=sig_char),stat="identity", color="gray20",size=0.1) +
  #   scale_fill_manual(values = color_vec_fullname, na.value="black") +
  #   theme_bw() +
  #   geom_hline(yintercept = 0) +
  #   scale_y_continuous(limits=c(-1,1)
  #                      # breaks = y_axis_values,
  #                      # labels = abs(y_axis_values)
  #   )+
  #   theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
  #   labs(title= plot_title,
  #        x="Average signature weight",
  #        # x="Variant",
  #        fill = "Signature") +
  #   theme(text = element_text(size=text_font_size),
  #         axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  #   guides(fill=guide_legend(ncol=legend_columns)) +
  #   theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  # 
  
  
  
  avg_weight_plot_out <- ggplot(avg_weights_sig_full_forplot) +
    geom_bar(aes(x=1,y = avg_weight,fill=sig_char),stat="identity", color="gray20",size=0.1) +
    coord_flip() +
    scale_fill_manual(values = color_vec_fullname, na.value="black") +
    theme_bw() +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits=c(-1,1),
                       # breaks = y_axis_values,
                       labels = c("1.0","0.5","0.0","0.5","1.0")
    )+
    theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=1)) +
    labs(title = plot_title,
         y="Average signature weight",
         # x="Variant",
         fill = "Signature") +
    theme(text = element_text(size=text_font_size),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
    guides(fill=guide_legend(ncol=legend_columns)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
  
  
  
  return(list(effect_plot = driver_plot_out + theme(legend.position = "none"),
              avg_weight_plot = avg_weight_plot_out+ theme(legend.position = "none")))
  
}









library(cowplot)




## Driver-centric figures 



# BLCA ---- 
blca_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "BLCA",
                                                      Bailey_tumor_type = "BLCA",
                                                      signatures_below_line = c("SBS2","SBS13"),
                                                      signatures_below_line_grouped = "APOBEC (2,13)",
                                                      drivers_to_plot = 10,
                                                      ordered_by_total_weight = F,
                                                      text_font_size = 10,
                                                      choose_variants_by_total_volume = T, 
                                                      plot_title = "BLCA")

plot_grid(blca_output$avg_weight_plot + labs(title=NULL), 
          blca_output$effect_plot + theme(text=element_text(size=10)),
          nrow = 1,rel_widths = c(1.1,5),align = "h")

blca_plot <- plot_grid(blca_output$avg_weight_plot + 
                         theme(axis.title.x = element_blank()) + 
                         theme(text=element_text(size=15)), 
                       blca_output$effect_plot + labs(title=NULL) + 
                         theme(text=element_text(size=15)) + 
                         theme(axis.title.x = element_blank()),
                       nrow = 2,rel_heights =  c(1.3,5),align = "v")


ggsave(filename = "dev/BLCA_driver_focus_total_volume_sort.png",height = 4,width = 5,plot = blca_plot)
# 
# 
# ```{r LUAD}
# 
luad_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "LUAD",
                                                      Bailey_tumor_type = "LUAD",
                                                      signatures_below_line = c("SBS4"),
                                                      signatures_below_line_grouped = c("Tobacco (4,29)"),
                                                      drivers_to_plot = 10,
                                                      ordered_by_total_weight = F,
                                                      text_font_size = 10,choose_variants_by_total_volume = T,
                                                      plot_title = "LUAD")
# 
# plot_grid(luad_output$avg_weight_plot + labs(title=NULL), luad_output$effect_plot,nrow = 1,rel_widths = c(1,5),align = "h")
# 

luad_plot <- plot_grid(luad_output$avg_weight_plot + 
                         theme(axis.title.x = element_blank()) + 
                         theme(text=element_text(size=15)), 
                       luad_output$effect_plot + labs(title=NULL) + 
                         theme(text=element_text(size=15)) + 
                         theme(axis.title.x = element_blank()),
                       nrow = 2,rel_heights =  c(1.7,5),align = "v")
ggsave(filename = "dev/LUAD_driver_focus_total_volume_sort.png",height = 3.5,width = 5,
       luad_plot)
# ```
# 
# 


# CESC ---- 

cesc_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "CESC",
                                                      Bailey_tumor_type = "CESC",
                                                      signatures_below_line = c("SBS2","SBS13"),
                                                      signatures_below_line_grouped = "APOBEC (2,13)",
                                                      drivers_to_plot = 10,
                                                      ordered_by_total_weight = F,
                                                      text_font_size = 10,
                                                      choose_variants_by_total_volume = T, 
                                                      plot_title = "CESC")
# 
# plot_grid(luad_output$avg_weight_plot + labs(title=NULL), luad_output$effect_plot,nrow = 1,rel_widths = c(1,5),align = "h")
# 

# cesc_plot <- plot_grid(cesc_output$avg_weight_plot + 
#                          theme(axis.title.x = element_blank()) + 
#                          theme(text=element_text(size=15)), 
#                        cesc_output$effect_plot + labs(title=NULL) + 
#                          theme(text=element_text(size=15)) + 
#                          theme(axis.title.x = element_blank()),
#                        nrow = 2,rel_heights =  c(1.3,5),align = "v")
# ggsave(filename = "dev/CESC_driver_focus_total_volume_sort.png",height = 3.5,width = 5,
#        cesc_plot)
# 
# 
# 
# 
# 
# 
# BLCA_CESC <- plot_grid(blca_plot,cesc_plot,nrow=1, align = "h")
# 
# 
# ggsave(filename = "dev/BLCA_CESC.png",height = 4,width = 8,
#        BLCA_CESC)





# BLCA_test <- plot_grid(blca_output$avg_weight_plot + 
#                          theme(axis.title.x = element_blank(),
#                                axis.text=element_text(size=12)) + 
#                          theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
#                        blca_output$effect_plot + 
#                          theme(axis.title.x = element_blank(), 
#                                title = element_blank(),
#                                axis.title.y = element_blank(),
#                                axis.text = element_text(size=12))+
#                          theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
#                        nrow = 2, align = "hv",
#                        rel_heights = c(1,3) )




BLCA_CESC <- plot_grid(blca_output$avg_weight_plot + 
                         theme(axis.title.x = element_blank(),
                               axis.text=element_text(size=12)) + 
                         theme(plot.margin = unit(c(0, 5, 0, 0), "pt")), 
                       cesc_output$avg_weight_plot + 
                         theme(axis.title.x = element_blank(),
                               axis.text=element_text(size=12)) + 
                         theme(plot.margin = unit(c(0, 5, 0, 0), "pt")), 
                       blca_output$effect_plot + 
                         theme(axis.title.x = element_blank(), 
                               title = element_blank(),
                               axis.title.y = element_blank(),
                               axis.text = element_text(size=12))+
                         theme(plot.margin = unit(c(0, 5, 0, 0), "pt")), 
                       cesc_output$effect_plot + 
                         theme(axis.title.x = element_blank(),
                               title = element_blank(),
                               axis.title.y = element_blank(),
                               axis.text=element_text(size=12)) + 
                         theme(plot.margin = unit(c(0, 5, 0, 0), "pt")),
                       nrow = 2, 
                       align = "v",rel_heights = c(1,3),labels = c("E","F") )


ggsave(filename = "dev/BLCA_CESC.png",height = 4,width = 8,
       BLCA_CESC)




# HPV + - -------


HNSC_pos_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "HNSC_HPVpos",
                                                          Bailey_tumor_type = "HNSC",
                                                          signatures_below_line = c("SBS2","SBS13"),
                                                          signatures_below_line_grouped = "APOBEC (2,13)",
                                                          drivers_to_plot = 10,
                                                          ordered_by_total_weight = F,
                                                          text_font_size = 10,
                                                          choose_variants_by_total_volume = T, 
                                                          plot_title = "HNSC HPV positive")



HNSC_neg_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "HNSC_HPVneg",
                                                          Bailey_tumor_type = "HNSC",
                                                          signatures_below_line = c("SBS2","SBS13"),
                                                          signatures_below_line_grouped = "APOBEC (2,13)",
                                                          drivers_to_plot = 10,
                                                          ordered_by_total_weight = F,
                                                          text_font_size = 10,
                                                          choose_variants_by_total_volume = T, 
                                                          plot_title = "HNSC HPV negative")




HNSC_combined <- plot_grid(HNSC_neg_output$avg_weight_plot + 
                             theme(axis.title.x = element_blank(),
                                   axis.text=element_text(size=12)) + 
                             theme(plot.margin = unit(c(5, 5, 0, 0), "pt")), 
                           HNSC_pos_output$avg_weight_plot + 
                             theme(axis.title.x = element_blank(),
                                   axis.text=element_text(size=12)) + 
                             theme(plot.margin = unit(c(5, 5, 0, 0), "pt")), 
                           HNSC_neg_output$effect_plot + 
                             theme(axis.title.x = element_blank(), 
                                   title = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.text = element_text(size=12))+
                             theme(plot.margin = unit(c(0, 5, 0, 0), "pt")), 
                           HNSC_pos_output$effect_plot + 
                             theme(axis.title.x = element_blank(),
                                   title = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.text=element_text(size=12)) + 
                             theme(plot.margin = unit(c(0, 5, 0, 0), "pt")),
                           nrow = 2, 
                           align = "v",rel_heights = c(1,3), labels = c("G","H") )


ggsave(filename = "dev/HNSC_combined.png",height = 4,width = 8,
       HNSC_combined)



# LUAD LUSC ------



LUAD_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "LUAD",
                                                      Bailey_tumor_type = "LUAD",
                                                      signatures_below_line = c("SBS4","SBS29"),
                                                      signatures_below_line_grouped = "Tobacco (4,29)",
                                                      drivers_to_plot = 10,
                                                      ordered_by_total_weight = F,
                                                      text_font_size = 10,
                                                      choose_variants_by_total_volume = T, 
                                                      plot_title = "LUAD")



LUSC_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "LUSC",
                                                      Bailey_tumor_type = "LUSC",
                                                      signatures_below_line = c("SBS4","SBS29"),
                                                      signatures_below_line_grouped = "Tobacco (4,29)",
                                                      drivers_to_plot = 10,
                                                      ordered_by_total_weight = F,
                                                      text_font_size = 10,
                                                      choose_variants_by_total_volume = T, 
                                                      plot_title = "LUSC")




Lung_combined <- plot_grid(LUAD_output$avg_weight_plot + 
                             theme(axis.title.x = element_blank(),
                                   axis.text=element_text(size=12)) + 
                             theme(plot.margin = unit(c(0, 5, 0, 0), "pt")), 
                           LUSC_output$avg_weight_plot + 
                             theme(axis.title.x = element_blank(),
                                   axis.text=element_text(size=12)) + 
                             theme(plot.margin = unit(c(0, 5, 0, 0), "pt")), 
                           LUAD_output$effect_plot + 
                             theme(axis.title.x = element_blank(), 
                                   title = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.text = element_text(size=12))+
                             theme(plot.margin = unit(c(0, 5, 0, 0), "pt")), 
                           LUSC_output$effect_plot + 
                             theme(axis.title.x = element_blank(),
                                   title = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.text=element_text(size=12)) + 
                             theme(plot.margin = unit(c(0, 5, 0, 0), "pt")),
                           nrow = 2, 
                           align = "v",rel_heights = c(1,3),labels = c("A","B") )


ggsave(filename = "dev/Lung_combined.png",height = 4,width = 8,
       Lung_combined)


# ESCA ---- 

ESCA_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "ESCA",
                                                      Bailey_tumor_type = "ESCA",
                                                      signatures_below_line = c("SBS16"),
                                                      signatures_below_line_grouped = "Alcohol-associated (16)",
                                                      drivers_to_plot = 10,
                                                      ordered_by_total_weight = F,
                                                      text_font_size = 10,
                                                      choose_variants_by_total_volume = T, 
                                                      plot_title = "ESCA")



# SKCM p m --------




SKCMp_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "SKCM_primary",
                                                       Bailey_tumor_type = "SKCM",
                                                       signatures_below_line = c("SBS7a",
                                                                                 "SBS7b",
                                                                                 "SBS7c",
                                                                                 "SBS7d",
                                                                                 "SBS38"),
                                                       signatures_below_line_grouped = "UV light (7a–d,38)",
                                                       drivers_to_plot = 10,
                                                       ordered_by_total_weight = F,
                                                       text_font_size = 10,
                                                       choose_variants_by_total_volume = T, 
                                                       plot_title = "SKCM (Primary)")



SKCMm_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "SKCM_metastasis",
                                                       Bailey_tumor_type = "SKCM",
                                                       signatures_below_line = c("SBS7a",
                                                                                 "SBS7b",
                                                                                 "SBS7c",
                                                                                 "SBS7d",
                                                                                 "SBS38"),
                                                       signatures_below_line_grouped = "UV light (7a–d,38)",
                                                       drivers_to_plot = 10,
                                                       ordered_by_total_weight = F,
                                                       text_font_size = 10,
                                                       choose_variants_by_total_volume = T, 
                                                       plot_title = "SKCM (Metastases)")




Skin_combined <- plot_grid(SKCMp_output$avg_weight_plot + 
                             theme(axis.title.x = element_blank(),
                                   axis.text=element_text(size=12)) + 
                             theme(plot.margin = unit(c(5, 5, 0, 0), "pt")), 
                           SKCMm_output$avg_weight_plot + 
                             theme(axis.title.x = element_blank(),
                                   axis.text=element_text(size=12)) + 
                             theme(plot.margin = unit(c(5, 5, 0, 0), "pt")), 
                           SKCMp_output$effect_plot + 
                             theme(axis.title.x = element_blank(), 
                                   title = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.text = element_text(size=12))+
                             theme(plot.margin = unit(c(0, 5, 0, 0), "pt")), 
                           SKCMm_output$effect_plot + 
                             theme(axis.title.x = element_blank(),
                                   title = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.text=element_text(size=12)) + 
                             theme(plot.margin = unit(c(0, 5, 0, 0), "pt")),
                           nrow = 2, 
                           align = "v",rel_heights = c(1,3),labels = c("G","H") )



ggsave(filename = "dev/Skin_combined.png",height = 4,width = 8,
       Skin_combined)




# BRCA ---- 

BRCA_ERp_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "BRCA_ER_pos",
                                                          Bailey_tumor_type = "BRCA",
                                                          signatures_below_line = c("SBS2","SBS13"),
                                                          signatures_below_line_grouped = "APOBEC (2,13)",
                                                          drivers_to_plot = 10,
                                                          ordered_by_total_weight = F,
                                                          text_font_size = 10,
                                                          choose_variants_by_total_volume = T, 
                                                          plot_title = "BRCA ER positive")





#  LIHC----
LIHC_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "LIHC",
                                                      Bailey_tumor_type = "LIHC",
                                                      signatures_below_line = c("SBS22",
                                                                                "SBS24",
                                                                                "SBS42",
                                                                                "SBS88"),
                                                      signatures_below_line_grouped = "Mutagenic chemical exposure (22,24,42,88)",
                                                      drivers_to_plot = 10,
                                                      ordered_by_total_weight = F,
                                                      text_font_size = 10,
                                                      choose_variants_by_total_volume = T, 
                                                      plot_title = "LIHC")


# OV ---- 

OV_output <- driver_focused_nrsi_weight_barplot_avg(tumor_type_ourdata = "OV",
                                                    Bailey_tumor_type = "OV",
                                                    signatures_below_line = c("SBS3"),
                                                    signatures_below_line_grouped = "Defective homologous recombination (3)",
                                                    drivers_to_plot = 10,
                                                    ordered_by_total_weight = F,
                                                    text_font_size = 10,
                                                    choose_variants_by_total_volume = T, 
                                                    plot_title = "OV")





skin_liver_combined <- plot_grid(SKCMp_output$avg_weight_plot + 
                                   theme(axis.title.x = element_blank(),
                                         axis.text=element_text(size=12)) + 
                                   theme(plot.margin = unit(c(5, 5, 0, 0), "pt")), 
                                 LIHC_output$avg_weight_plot + 
                                   theme(axis.title.x = element_blank(),
                                         axis.text=element_text(size=12)) + 
                                   theme(plot.margin = unit(c(5, 5, 0, 0), "pt")), 
                                 SKCMp_output$effect_plot + 
                                   theme(axis.title.x = element_blank(), 
                                         title = element_blank(),
                                         axis.title.y = element_blank(),
                                         axis.text = element_text(size=12))+
                                   theme(plot.margin = unit(c(0, 5, 0, 0), "pt")), 
                                 LIHC_output$effect_plot + 
                                   theme(axis.title.x = element_blank(),
                                         title = element_blank(),
                                         axis.title.y = element_blank(),
                                         axis.text=element_text(size=12)) + 
                                   theme(plot.margin = unit(c(0, 5, 0, 0), "pt")),
                                 nrow = 2, 
                                 align = "v",rel_heights = c(1,3),labels = c("C","D") )





all_driver_figure <- plot_grid(BLCA_CESC,
                               HNSC_combined,
                               Lung_combined,
                               skin_liver_combined,align = "v",ncol = 1)



cowplot::save_plot(all_driver_figure,filename = "dev/combined_driver_figure.png",base_height = 12,base_width = 8)





# legend work ---- 


leg_df <- tibble(x=1, Signature = factor(names(color_vec_fullname)))
leg_df$Signature <- fct_relevel(leg_df$Signature, names(color_vec_fullname))

leg_ggplot <- ggplot(leg_df, aes(x=x,y=1,fill=Signature)) + 
  geom_col() + 
  scale_fill_manual(values = color_vec_fullname) + 
  theme(legend.position="bottom") + 
  guides(fill=guide_legend(ncol=2))

leg_legend <- suppressWarnings( cowplot::get_legend(leg_ggplot))



all_driver_figure_w_legend <- plot_grid(all_driver_figure,leg_legend, ncol=1,rel_heights = c(1,0.1))


cowplot::save_plot(all_driver_figure_w_legend,filename = "dev/combined_driver_figure_w_legend.png",base_height = 14,base_width = 8)



# all driver 2x2-----

all_driver_figure_2x2 <- plot_grid(
  Lung_combined,
  skin_liver_combined,
  BLCA_CESC,
  HNSC_combined,align = "v",ncol = 2)



cowplot::save_plot(all_driver_figure_2x2,filename = "dev/combined_driver_figure2x2.png",base_height = 7,base_width = 13)





# legend work2 ---- 


leg_df <- tibble(x=1, Signature = factor(names(color_vec_fullname)))
leg_df$Signature <- fct_relevel(leg_df$Signature, names(color_vec_fullname))

leg_ggplot <- ggplot(leg_df, aes(x=x,y=1,fill=Signature)) + 
  geom_col( color="gray20",size=0.1) + 
  scale_fill_manual(values = color_vec_fullname) + 
  theme(legend.position="bottom") + 
  guides(fill=guide_legend(nrow=2))

leg_legend <- suppressWarnings( cowplot::get_legend(leg_ggplot))




y_label <- ggdraw() +
  draw_label(
    "Variant",
    fontface = 'bold',
    x = 0,
    angle=90,
    y=0.5,
    vjust = 0.5,
    hjust = 0.5,
    size = 15
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

x_label <- ggdraw() +
  draw_label(
    "Proportion of mutational weight or cancer effect attributable to each mutational process",
    fontface = 'bold',
    x = 0.5175,
    angle=0,
    y=0,
    vjust = 1,
    hjust = 0.5,
    size = 15
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(7, 0, 20, 0)
  )



driver_fig_w_lab <- plot_grid(
  y_label, all_driver_figure_2x2,
  ncol = 2,
  # rel_heights values control vertical title margins
  rel_widths = c(0.015, 1)
)


driver_fig_w_lab <- plot_grid(
  driver_fig_w_lab, x_label,
  nrow = 2,
  # rel_heights values control vertical title margins
  rel_heights = c(1, 0.035)
)


all_driver_figure_w_legend_2x2 <- plot_grid(leg_legend,driver_fig_w_lab, ncol=1,rel_heights = c(0.1,1))







cowplot::save_plot(all_driver_figure_w_legend_2x2,filename = "dev/combined_driver_figure_w_legend_2x2.png",base_height = 7.5,base_width = 13)

cowplot::save_plot(all_driver_figure_w_legend_2x2,filename = "output_data/figure_3_variant_effects.png",base_height = 7.5,base_width = 13)


