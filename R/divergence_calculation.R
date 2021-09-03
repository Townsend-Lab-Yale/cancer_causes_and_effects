
library(tidyverse)
library(philentropy)

divergence_calculation <- function(trinuc_weights, NRSI_per_tumor){
  
  NRSI_data <- NRSI_per_tumor %>% 
    dplyr::select(variant,Unique_Patient_Identifier,ends_with("nrsi")) %>% # select nrsi values
    pivot_longer(cols = ends_with("nrsi")) %>% # make long
    group_by(Unique_Patient_Identifier) %>% 
    mutate(pct = value / sum(value)) %>% # nrsi is a percent of the nrsi for that tumor
    group_by(Unique_Patient_Identifier, name) %>%
    summarize(Weight_nrsi = sum(pct)) %>%
    mutate(Tumor = Unique_Patient_Identifier, 
           Signature = gsub(pattern = "_nrsi",replacement = "",x = name)) %>%
    ungroup() %>%
    dplyr::select(Tumor, Signature, Weight_nrsi)
  

  signature_weights <- trinuc_weights %>% 
    dplyr::select(-total_snvs,-sig_extraction_snvs,-group_avg_blended) %>%
    pivot_longer(-Unique_Patient_Identifier) %>%
    mutate(Tumor = Unique_Patient_Identifier, 
           Signature = name, 
           Weight = value) %>%
    dplyr::select(-Unique_Patient_Identifier, -name, -value) %>%
    group_by(Tumor) %>%
    mutate(Weight_prop = Weight / (sum(Weight))) %>%
    ungroup()
  
  
  
  NRSI_data_complete <- NRSI_data %>%
    tidyr::complete(Signature = unique(signature_weights$Signature),
                    Tumor,
                    fill=list(Weight_nrsi = 0))
  
  
  
  # start collapse block -----
  
  
  signatures_names_matrix <- cancereffectsizeR::get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.1")
  signatures_names_matrix <- as.data.frame(signatures_names_matrix$meta[,c("Signature","Etiology")])
  
  signatures_names_matrix$Etiology_sig <- paste0(signatures_names_matrix$Etiology, " (",gsub(pattern = "SBS",replacement = "",x = signatures_names_matrix$Signature),")")
  
  rownames(signatures_names_matrix) <- signatures_names_matrix$Signature
  
  
  signatures_names_matrix$sig_num <- gsub(pattern = "SBS",replacement = "",x = signatures_names_matrix$Signature)
  
  
  NRSI_data_complete %>%
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
    NRSI_data_complete
  
  
  
  signature_weights <- signature_weights %>% mutate(signature_full = signatures_names_matrix[as.character(Signature), "Etiology_sig"]) %>%
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
                                                  `Alcohol-associated (16)` = "Unknown (16)")))
  
  
  signature_weights <- signature_weights %>% 
    group_by(Tumor,signature_full) %>% 
    summarize(Weight_prop = sum(Weight_prop))
  
  NRSI_data_complete <- NRSI_data_complete %>% 
    group_by(Tumor,signature_full) %>% 
    summarize(Weight_nrsi = sum(Weight_nrsi))
  
  
  # end collapse block -----
  
  tumor_JSD_df <- data.frame(Tumor = as.character(unique(NRSI_data_complete$Tumor)),JSD=NA,stringsAsFactors = F)
  
  for(tumor_ind in 1:length(unique(NRSI_data_complete$Tumor))){
    
    NRSIs <- NRSI_data_complete[as.character(NRSI_data_complete$Tumor) == tumor_JSD_df$Tumor[tumor_ind],]
    sig_weights <- signature_weights[signature_weights$Tumor == tumor_JSD_df$Tumor[tumor_ind],]
    
    merged_NRSIs_sig_weights <- merge(NRSIs, sig_weights)
    
    tumor_JSD_df$JSD[tumor_ind] <- philentropy::JSD(x = rbind(merged_NRSIs_sig_weights$Weight_nrsi,
                                                              merged_NRSIs_sig_weights$Weight_prop))
    
  }
  
  
  return(list(tumor_JSD_df=tumor_JSD_df,
              tumor_signature_weights=signature_weights))
  
}


