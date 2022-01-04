

# fig4 data gather
library(tidyverse)
library(data.table)


file_names <- paste0("combined_attribution_")

file_names <- dir(recursive = T,pattern = paste0(file_names,collapse = "|"),full.names = T)

# UCEC takes up too much memory, needs stand alone solution
file_names <- file_names[str_detect(string = file_names,pattern = "UCEC")]

## for dev
# file_names <- "./bootstrap_analysis/combined_attribution_LUAD.rds"


# objects to hold data summaries
fig_4_storage <- vector(mode = "list",length = length(file_names))
names(fig_4_storage) <- file_names


for(file_ind in seq_along(file_names)){
  
  # file_ind <- 1
  # effect size attribution----- 
  these_attributions <- readRDS(file = file_names[file_ind])
  these_attributions <- rbindlist(these_attributions, idcol = "bootstrap_sample",fill = T)
  
  message("Loaded in the data")
  
  # different tumors had different attribution columns, fill in NA above then replace with zero (there was zero attribution to that missing signature if missing)
  these_attributions <- these_attributions %>%
    replace(is.na(.),0)
  
  message("Added in the zeroes")
  
  tumors_w_attribution <- unique(these_attributions$Unique_Patient_Identifier)
  
  this_file <- file_names[file_ind]
  subnames <- str_split(this_file,pattern = ".rds")[[1]]
  subnames <- str_split(subnames[1],pattern = "combined_attribution_")[[1]]
  tumor_type <- subnames[length(subnames)]
  
  
  # effect_size_attribution <- these_attributions %>%
  #   # do not need variant information because we are summing across variants to find signature contribution
  #   select(-ends_with("flux"),-variant) %>%
  #   pivot_longer(cols = starts_with("SBS"), names_to = "signature",values_to = "weight") %>%
  #   group_by(bootstrap_sample, Unique_Patient_Identifier) %>%
  #   mutate(proportional_weight = weight / sum(weight)) %>%
  #   ungroup() %>%
  #   group_by(bootstrap_sample, signature) %>%
  #   summarize(total_effect_contribution = sum(proportional_weight)/length(tumors_w_attribution)) %>%
  #   group_by(signature) %>%
  #   summarize(avg_total_effect_contribution = mean(total_effect_contribution))
  #   
  # 
  # effect_size_attribution <- effect_size_attribution %>%
  #   mutate(signature = stringr::str_remove(string = signature, 
  #                                          pattern = "_nrsi")) %>% 
  #   mutate(data_type = "effect_size_attribution") %>%
  #   mutate(tumor_type = tumor_type)
  
  
  these_attributions <- these_attributions %>%
    select(-ends_with("flux")) %>%
    pivot_longer(cols = starts_with("SBS"), names_to = "signature",values_to = "weight")
  
  message("pivoted first")
  gc()
  
  these_attributions <- these_attributions %>%
    group_by(bootstrap_sample, Unique_Patient_Identifier) %>%
    mutate(proportional_weight = weight / sum(weight)) %>%
    select(-weight) %>%
    ungroup() 
  
  message("first mutate")
    
  these_attributions <- these_attributions %>%
    group_by(bootstrap_sample, Unique_Patient_Identifier,signature) %>% 
    summarize(proportional_weight = sum(proportional_weight)) %>%
    group_by(Unique_Patient_Identifier, signature) %>%
    summarize(avg_proportional_weight = mean(proportional_weight)) %>%
    ungroup()
  
  gc()
  
  message("finished attribution analyses")
  
  # trinuc_attribution ----- 
  
  trinuc_file <- stringr::str_replace(string = file_names[file_ind],
                                      pattern = "combined_attribution_",
                                      replacement = "combined_trinuc_storage_")
  
  these_trinuc <- readRDS(file = trinuc_file)
  these_trinuc <- rbindlist(these_trinuc,idcol="bootstrap_sample")
  
  # trinuc_avg_weight <- these_trinuc %>% 
  #   filter(Unique_Patient_Identifier %in% tumors_w_attribution) %>% #only want tumors where we have attribution data, i.e. recurrent variants
  #   select(-ends_with("snvs"),-group_avg_blended) %>%
  #   pivot_longer(cols = starts_with("SBS"), names_to = "signature",values_to = "weight") %>% 
  #   filter(weight > 0) %>%
  #   group_by(bootstrap_sample,signature) %>% 
  #   summarize(total_weight = sum(weight)/length(tumors_w_attribution)) %>% 
  #   group_by(signature) %>% 
  #   summarize(avg_total_weight = mean(total_weight))
  # 
  
  these_trinuc <- these_trinuc %>% 
    filter(Unique_Patient_Identifier %in% tumors_w_attribution) %>% #only want tumors where we have attribution data, i.e. recurrent variants
    select(-ends_with("snvs"),-group_avg_blended) %>%
    pivot_longer(cols = starts_with("SBS"), names_to = "signature",values_to = "weight") %>%
    group_by(Unique_Patient_Identifier, signature) %>% 
    summarize(avg_proportional_weight = mean(weight)) %>% ungroup()
  
  gc()
  # combine ----- 
  
  these_attributions <- these_attributions %>% 
    mutate(signature = stringr::str_remove(string = signature, pattern = "_nrsi"))
  
  joined_df <- left_join(x = these_trinuc,
                         y = these_attributions,
                         by = c("Unique_Patient_Identifier","signature"),
                         suffix = c(".trinuc",".effectsize") ) 
  
  joined_df %>%
    mutate(avg_proportional_weight.effectsize = case_when(
      is.na(avg_proportional_weight.effectsize) ~ 0,
      TRUE ~ avg_proportional_weight.effectsize)) %>%
    pivot_longer(cols = starts_with("avg"),values_to = "weight",names_to = "data_type") %>% 
    mutate(data_type = stringr::str_remove(string = data_type, pattern = "avg_proportional_weight.")) -> 
    joined_df
  
  
  
  fig_4_storage[[file_ind]] <- joined_df
  names(fig_4_storage)[file_ind] <- tumor_type
  
  rm(these_attributions,these_trinuc)
  gc()
  
  print(file_ind) 
  
}

# ran with LUAD to test UCEC before splitting by list 
# saveRDS(object = fig_4_storage,file = "dev/fig_4_storage_LUAD_ucectest.rds")

saveRDS(object = fig_4_storage,file = "fig_4_storage_UCEC.rds")





