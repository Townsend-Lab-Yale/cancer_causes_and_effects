

# define color palette for the obects in manuscript with colors 
# thanks https://colorbrewer2.org/ 


test_colors <- tribble(
  ~ object,                               ~ objects_color, 
  "KLF5 E419Q",                           "#377eb8",
  "OR2T34 L163L",                         "#4daf4a",
  "TP53 R282W",                           "#e41a1c",
  "Deamination with age, clock-like (1)", "gray40",
  "Unknown, clock-like (5)",              "gray60",
  "APOBEC (2,13)",                        "#7570b3",
  "Defective homologous recombination (3)","#e7298a",
  "Tobacco (4,29)",                       "#a6761d",
  "UV light (7a–d,38)",                   "#e6ab02",
  "Prior treatment (11,31,32,35)",         "#1b9e77",
  "Mutagenic chemical exposure (22,24,42,88)", "#66a61e",
  "Alcohol-associated (16)",                "#d95f02",
  "Non-actionable and unknown signatures", "black"
  
)


test_colors$object <- factor(test_colors$object,levels=rev(test_colors$object))

color_vec <- setNames(nm= as.character(test_colors$object),object = as.character(test_colors$objects_color))


# Testing how it looks ----
# ggplot(test_colors, aes(x=object,y=1,fill=object)) +
#   geom_col() +
#   scale_fill_manual(values = color_vec) +
#   theme_classic() +
#   coord_flip()
# 
# 
# ggsave(filename = "manuscript_analysis_bootstrap/scratch/color_palette_example.png")
# # 
# 
# 

# converting to SBS for plots 

signatures <- cancereffectsizeR::get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.2")
signatures <- signatures$meta$Signature

sig_tib <- tibble(signature = signatures) 


sig_tib <- sig_tib  %>%
  mutate(color_choice = case_when(
    signature == "SBS1" ~ as.character(color_vec["Deamination with age, clock-like (1)"]),
    signature == "SBS5" ~ as.character(color_vec["Unknown, clock-like (5)"]),
    signature %in% c("SBS2","SBS13") ~ as.character(color_vec["APOBEC (2,13)"]),
    signature %in% c("SBS3") ~ as.character(color_vec["Defective homologous recombination (3)"]),
    signature %in% c("SBS4","SBS29") ~ as.character(color_vec["Tobacco (4,29)"]),
    signature %in% c("SBS7a","SBS7b","SBS7c","SBS7d","SBS38") ~ as.character(color_vec["UV light (7a–d,38)"]),
    signature %in% c("SBS11","SBS31","SBS32","SBS35") ~ as.character(color_vec["Prior treatment (11,31,32,35)"]),
    signature %in% c("SBS22","SBS24","SBS42","SBS88") ~ as.character(color_vec["Mutagenic chemical exposure (22,24,42,88)"]),
    signature %in% c("SBS16") ~ as.character(color_vec["Alcohol-associated (16)"]),
    TRUE ~ color_vec["Non-actionable and unknown signatures"]))




color_vec_sbs <- setNames(nm= as.character(sig_tib$signature),object = as.character(sig_tib$color_choice))
