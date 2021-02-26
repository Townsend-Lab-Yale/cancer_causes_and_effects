# script to combine all effect size and population-scaled effect size data

library(tidyverse)
library(cancereffectsizeR)


# combining effect size---- 
message("Combining effect size data...")


effect_size_files <-  dir(full.names = T,recursive = T)[
  grep(pattern = "__selection_output.RData",x = dir(recursive = T))]



tumor_names <- unlist(strsplit(effect_size_files,split = "/"))[c(F,T,F)]



selection_data_main <- NULL
variant_prevalence_main <- NULL
total_substitutions <- vector(mode = "list",length = length(tumor_names))
names(total_substitutions) <- tumor_names

for(j in 1:length(effect_size_files)){
  
  
  load(effect_size_files[j])
  
  
  maf <- analysis@maf[variant_type=="snv"]
  total_substitutions[[j]] <- table(maf$Unique_Patient_Identifier)
  
  if("total_maf_freq" %in% colnames(analysis$variants)){
    variant_prevalence  <- as.data.frame(analysis$variants) %>%
      mutate(maf_frequency = total_maf_freq) %>%
      mutate(samples_covering = total_maf_freq) %>%
      filter(maf_frequency > 1) %>%
      select(variant_name,samples_covering,maf_frequency) %>%
      mutate(tumor_type = tumor_names[j])
  }else{
  
  variant_prevalence <- analysis$variants %>%
    filter(maf_frequency > 1) %>%
    select(variant_name,samples_covering,maf_frequency) %>%
    mutate(tumor_type = tumor_names[j])
  
  }
  
  variant_prevalence_main <- rbind(variant_prevalence_main, variant_prevalence)
  
  selection_data <- as.data.frame(analysis@selection_results$selection.1)
  
  
  selection_data$tumor_type <- tumor_names[j]
  
  
  
  selection_data_main <- rbind(selection_data_main,selection_data)
  
  
  # View(selection_data_df)
  print(tumor_names[j])
  
}

save(selection_data_main,file = "combined_selection_results.RData")

save(total_substitutions, file = "combined_substitution_data.RData")

save(variant_prevalence_main, file = "variant_prevalence_main.RData")


# combining JSD ---- 
message("Combining JSD data...")


JSD_files <- dir(full.names = T,recursive = T)[grep(pattern = "_tumor_JSD_df.RData",x = dir(recursive = T))]



tumor_names <- unlist(strsplit(JSD_files,split = "/"))[c(F,T,F)]

JSD_combined <- NULL
weights_combined <- NULL
for(ind in 1:length(JSD_files)){
  
  load(JSD_files[ind])
  
  tumor_JSD_df$tumor_JSD_df$tumor_type <- tumor_names[ind]
  tumor_JSD_df$tumor_signature_weights$tumor_type <- tumor_names[ind]
  
  JSD_combined <- rbind(JSD_combined,tumor_JSD_df$tumor_JSD_df)
  weights_combined <- rbind(weights_combined,tumor_JSD_df$tumor_signature_weights)
  
  print(tumor_names[ind])
}

save(JSD_combined, file="JSD_combined.RData")
save(weights_combined, file="weights_combined.RData")


# per tumor scaled selection and combine -----
message("Calculating the per tumor scaled selection and combining data...")

per_tumor_data_files <- dir(full.names = T,recursive = T)[grep(pattern = "_population_scaled_selection.RData",x = dir(recursive = T))]

tumor_names <- unlist(strsplit(per_tumor_data_files,split = "/"))[c(F,T,F)]



per_tumor_data_main <- NULL
recurrent_var_per_tumor_main <- NULL

all_scaled_selection_main <- NULL

library(plyr)


for(j in 1:length(per_tumor_data_files)){
  
  
  load(per_tumor_data_files[j])
  
  all_scaled_selection_main <- plyr::rbind.fill(all_scaled_selection_main, 
                         population_scaled_selection %>%
                           mutate(tumor_type = tumor_names[j]))
  
  per_tumor_data <- population_scaled_selection %>%
    dplyr::select(variant,Unique_Patient_Identifier,ends_with("nrsi")) %>% # select nrsi values
    pivot_longer(cols = ends_with("nrsi"),
                 names_to = "Signature",values_to="NRSI") %>% # make long
    group_by(Unique_Patient_Identifier) %>% 
    mutate(NRSI_pct = NRSI / sum(NRSI))
  
  per_tumor_data$tumor_type <- tumor_names[j]
  
  
  
  per_tumor_data_main <- rbind(per_tumor_data_main,per_tumor_data)
  
  
  recurrent_var_per_tumor <- per_tumor_data %>%
    group_by(Unique_Patient_Identifier) %>%
    summarize(count=n_distinct(variant))
  
  
  recurrent_var_per_tumor$tumor_type <- tumor_names[j]
  
  recurrent_var_per_tumor_main <- rbind(recurrent_var_per_tumor_main,recurrent_var_per_tumor)
  
  # View(per_tumor_data_df)
  print(tumor_names[j])
}


save(per_tumor_data_main,file = "combined_per_tumor_proportion_scaled_selection.RData")

save(recurrent_var_per_tumor_main,file = "recurrent_var_per_tumor.RData")

save(all_scaled_selection_main, file="all_scaled_selection_main.RData")


