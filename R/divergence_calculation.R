
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


