
# figure showing trinucleotide context of signatures 


# just get the signatures with > 0 median weight
weight_summary %>% 
  filter(median>0) %>% 
  pull(signature) %>%
  as.character() -> 
  signatures_to_plot


# get our trinuc contexts 
# commented out because only need to run it once to find the context in these tumors. 

# # load in one of the cesa objects from the LUSC bootstrap run
# cesa_from_bootstrap <- cancereffectsizeR::load_cesa("bootstrap_analysis/bootstrap_results/LUSC/LUSC_cesa_1.rds")
# 
# # Coding mentorship from J. Mandell -> 
# # Let's say you have a vector of coding variants that you want trinuc_muts for:
# my_aac_muts = c("KLF5_E419Q", "OR2T34_L163L", "TP53_R282W")
# # calling select_variants with include_subvariants = T causes SNV records to be returned even when they are part of coding mutations
# # supplying variant_passlist filters results; you could leave out the option to get all annotations on all variant types
# variants = select_variants(cesa = cesa_from_bootstrap, variant_ids = my_aac_muts, include_subvariants = T)
# # subset to just the SNV records and you should see trinuc_mut for each, and the all_aac column identifies which coding variant each is associated with.
# snvs = variants[variant_type == "snv"]
# 
# snvs %>% unnest(covered_in)
# 
# cesa_from_bootstrap$maf %>% 
#   filter(Unique_Patient_Identifier == "TCGA-98-A53J-01A-11D-A26M-08") -> 
#   this_tumor_maf
# 
# snvs[snvs$variant_name %in% this_tumor_maf$variant_id,] %>% 
#   select(gene, trinuc_mut) -> 
#   variant_trinuc_contexts




# gene trinuc_mut
# 1: OR2T34    G[C>T]A
# 2:   KLF5    T[C>G]A
# 3:   TP53    C[C>T]G

sig_context <- cancereffectsizeR::get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.2")
sig_context <- sig_context$signatures[signatures_to_plot,]

# just want contexts relevant to our variants for the figure
sig_context <- sig_context[,stringr::str_detect(string = colnames(sig_context),pattern = "C>T|C>G")]

sig_context_long <- sig_context %>% 
  mutate(sig_names = rownames(sig_context)) %>% 
  pivot_longer(-sig_names,names_to = "trinuc",values_to = "weight")

sig_context_long <- sig_context_long %>%
  mutate(facet_labels = case_when(
    sig_names == "SBS1" ~ "Deamination with age (1)", 
    sig_names == "SBS2" ~ "APOBEC (2)", 
    sig_names == "SBS5" ~ "Unknown, clock-like (5)",
    sig_names == "SBS4" ~ "Tobacco (4)"
    ))



sig_context_long$facet_labels <- factor(sig_context_long$facet_labels,
                                        levels = unique(sig_context_long$facet_labels))

trinuc_colors <- setNames(object = rep("black",length(sig_context_long$trinuc)),nm = sig_context_long$trinuc)

# gene trinuc_mut
# 1: OR2T34    G[C>T]A
# 2:   KLF5    T[C>G]A
# 3:   TP53    C[C>T]G



trinuc_colors["G[C>T]A"] <- color_vec["OR2T34 L163L"]
trinuc_colors["T[C>G]A"] <- color_vec["KLF5 E419Q"]
trinuc_colors["C[C>T]G"] <- color_vec["TP53 R282W"]


sig_context_long$trinuc <- factor(sig_context_long$trinuc, levels = unique(sig_context_long$trinuc))

color_vec_axis_label <- tibble(trinuc=levels(sig_context_long$trinuc))

color_vec_axis_label$color <- trinuc_colors[color_vec_axis_label$trinuc]

sig_context_long %>% 
  ggplot(aes(x=trinuc,y=weight,fill=sig_names)) + 
  geom_col(aes(color=trinuc)) +
  scale_color_manual(values = trinuc_colors) + 
  facet_wrap(~facet_labels,nrow=4,scales = "free_y") + 
  theme_classic() + 
  scale_fill_manual(values = color_vec_sbs,limits=force) + 
  guides(fill="none",color="none") + 
  theme(axis.text.x = element_text(angle=90,vjust=.5)) + 
  labs(y="Trinucleotide weight", x="Trinucleotide mutation") +   
  theme(axis.text.x = element_text(color=color_vec_axis_label$color),
        text = element_text(size = plot_text_size)) -> 
  fig1_trinuc_definitions





