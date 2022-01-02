#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



# args <- "LUSC" # for debugging 

library(tidyverse)

source("../parameters_for_tumor_types.R")

parameters_for_run <- parameters_for_run[[ args[1] ]]

n_cores <- parameters_for_run %>% 
  filter(input_parameter == "n_cores") %>% 
  pull(input_value)

tumor_name <- parameters_for_run %>% 
  filter(input_parameter == "tumor_name") %>% 
  pull(input_value)

# combine and save trinuc and attribution files ---- 

# find the appropriate files
trinuc_files <- dir(path = "bootstrap_analysis",pattern = "trinuc_storage_[0-9]",recursive = T,full.names = T)
attribution_files <- dir(path = "bootstrap_analysis",pattern = "attribution_storage_[0-9]",recursive = T,full.names = T)

trinuc_storage <- vector(mode = "list",length = 1000)
attribution_storage <- vector(mode = "list", length = 1000)

# loop through files to include them all in a big list
for(file_ind in seq_along(trinuc_files)){
  
  this_trinuc <- readRDS(file = trinuc_files[file_ind])
  this_attribution <- readRDS(file = attribution_files[file_ind])
  
  filled_in <- which(lapply(this_trinuc,length)>0)
  
  trinuc_storage[filled_in] <- this_trinuc[filled_in]
  attribution_storage[filled_in] <- this_attribution[filled_in]
  
  
}

# save the list
saveRDS(object = trinuc_storage, file = paste0("bootstrap_analysis/combined_trinuc_storage_",tumor_name,".rds"))
saveRDS(object = attribution_storage, file = paste0("bootstrap_analysis/combined_attribution_",tumor_name,".rds"))



# saving selection data for figure 1 ------ 

if(tumor_name == "LUSC"){
  
  selection_storage <- vector(mode = "list",length = 1000)
  
  my_aac_muts <- c("KLF5_E419Q", "OR2T34_L163L", "TP53_R282W")
  cesa_files <- dir(path = "bootstrap_analysis/",pattern = "LUSC_cesa_[0-9]",recursive = T,full.names = T)
  
  
  for(cesa_ind in seq_along(cesa_files)){
    
    this_cesa <- cancereffectsizeR::load_cesa(file = cesa_files[cesa_ind])
    
    these_names <- names(this_cesa$selection)
    these_selections <- stringr::str_remove(string = names(this_cesa$selection),pattern = "boot_") %>% 
      as.numeric()
    
    names_and_selections <- data.frame(names = these_names,selections=these_selections)
    
    
    for(selection_ind in 1:nrow(names_and_selections)){
      
      this_cesa$selection[[names_and_selections[selection_ind,"names"]]] %>%
        filter(stringr::str_detect(string = variant_id,pattern = paste(my_aac_muts,collapse = "|"))) -> 
        selection_storage[[names_and_selections[selection_ind,"selections"]]]
      
      
    }
  }
  
  
  saveRDS(object = selection_storage, file = paste0("bootstrap_analysis/combined_selection_",tumor_name,".rds"))
  
}





# JSD calculations ---- 

library(data.table)
library(philentropy)

## examples from philentropy package 
# P <- 1:10/sum(1:10)
# Q <- 20:29/sum(20:29)
# x <- rbind(P,Q)
# suppressMessages(JSD(x))
# 
# 
# Prob <- rbind(1:10/sum(1:10), 20:29/sum(20:29), 30:39/sum(30:39))
# 
# # compute the KL matrix of a given probability matrix
# JSDMatrix <- JSD(Prob)


# need to perform JSD calculations on tumor-by-tumor basis
# recombining data structures so they are aligned by tumors and not bootstraps
trinuc_storage_singleDT <- rbindlist(trinuc_storage,idcol = "bootstrap_sample",fill = T)

trinuc_storage_singleDT <- trinuc_storage_singleDT %>% 
  replace(is.na(.),0)

trinuc_storage_list_by_tumor <- split(trinuc_storage_singleDT,f = trinuc_storage_singleDT$Unique_Patient_Identifier)

attribution_storage_singleDT <- rbindlist(attribution_storage, idcol = "bootstrap_sample",fill = T)

# different tumors had different attribution columns, fill in NA above then replace with zero (there was zero attribution to that missing signature if missing)
attribution_storage_singleDT <- attribution_storage_singleDT %>%
  replace(is.na(.),0)

attribution_storage_list_by_tumor <- split(attribution_storage_singleDT, f = attribution_storage_singleDT$Unique_Patient_Identifier)

rm(trinuc_storage_singleDT,attribution_storage_singleDT) # try and clean up some big data 
gc()


# will use parallel processing so need a big list to run through 
jsd_storage_list <- vector(mode = "list",length = length(attribution_storage_list_by_tumor))
names(jsd_storage_list) <- names(attribution_storage_list_by_tumor)

for(tumor_ind in seq_along(jsd_storage_list)){
  
  jsd_storage_list[[tumor_ind]] <- list(trinuc = trinuc_storage_list_by_tumor[[names(jsd_storage_list)[tumor_ind]]],
                                        attribut = attribution_storage_list_by_tumor[[names(jsd_storage_list)[tumor_ind]]])
  
  
}



# this_tumor <- jsd_storage_list$`MDA-1006-T` # debug 
options(dplyr.summarise.inform = FALSE)

jsd_per_tumor <- function(this_tumor){
  
  tumor_name <- this_tumor$trinuc$Unique_Patient_Identifier[1]
  
  # trinuc_long <- this_tumor$trinuc %>%
  #   select(-Unique_Patient_Identifier,-ends_with("_snvs"),-group_avg_blended) %>%
  #   pivot_longer(cols = starts_with("SBS"),names_to = "signature",values_to = "weight")
  # 
  
  
  # need to create two matrices with the values to pull from JSD
  # need them to have the same columns referencing the same signatures
  # and, need them to sum to 1. 
  trinuc_wide <- this_tumor$trinuc %>% select(starts_with("SBS"))
  
  # need to collapse signatures that are the same "mechanism" as defined in the manuscript
  trinuc_wide %>% 
    mutate(bootstrap_sample = 1:nrow(.)) %>%
    pivot_longer(cols = -bootstrap_sample, names_to = "signature",values_to = "weight") %>%
    mutate(signature_merged = as.factor(signature)) %>%
    mutate(signature_merged = forcats::fct_collapse(
      signature_merged, 
      `SBS7` = c("SBS7a","SBS7b","SBS7c","SBS7d"),
      `SBS2_13` = c("SBS2","SBS13"),
      `SBS4_29` = c("SBS4","SBS29"),
      `SBStreatment` = c("SBS11","SBS31","SBS32","SBS35"),
      `SBSmut` = c("SBS22","SBS24","SBS42","SBS88"))) %>%
    select(-signature) %>%
    group_by(bootstrap_sample,signature_merged) %>%
    summarize(weight = sum(weight)) %>%
    ungroup() %>%
    pivot_wider(names_from = signature_merged, values_from = weight) %>%
    select(-bootstrap_sample) -> 
    trinuc_wide
  
  
  
  ##  already sums to 1
  # all(!(rowSums(trinuc_wide)-1)>1e-15) # all very close to 1, within 1e-15
  
  attribute_wide <- this_tumor$attribut %>% 
    select(bootstrap_sample,ends_with("nrsi"))
  
  colnames(attribute_wide) <- stringr::str_remove(string = colnames(attribute_wide), pattern = "_nrsi")
  attribute_long <- attribute_wide %>% 
    pivot_longer(cols = -bootstrap_sample, names_to = "signature",values_to = "weight") %>%
    mutate(signature_merged = as.factor(signature)) %>%
    mutate(signature_merged = suppressWarnings(
      forcats::fct_collapse(
        signature_merged, 
        `SBS7` = c("SBS7a","SBS7b","SBS7c","SBS7d"),
        `SBS2_13` = c("SBS2","SBS13"),
        `SBS4_29` = c("SBS4","SBS29"),
        `SBStreatment` = c("SBS11","SBS31","SBS32","SBS35"),
        `SBSmut` = c("SBS22","SBS24","SBS42","SBS88")))) %>%
    mutate(signature = signature_merged) %>%
    group_by(bootstrap_sample,signature) %>%
    summarize(total_weight = sum(weight)) %>%
    group_by(bootstrap_sample) %>%
    mutate(proportion_weight = total_weight / sum(total_weight)) %>%
    ungroup()
  
  attribute_wide <- attribute_long %>% 
    select(-total_weight) %>%
    pivot_wider(names_from = signature,values_from =  proportion_weight) %>%
    select(-bootstrap_sample)
  
  # all(!(rowSums(attribute_wide)-1)>1e-15) # all very close to 1, within 1e-15
  
  attribute_long %>%
    pull(signature) %>%
    unique() -> 
    SBS_in_nrsi
  
  trinuc_wide <- trinuc_wide %>%
    select(all_of(SBS_in_nrsi))
  
 
  
  trinuc_wide <- as.matrix(trinuc_wide)
  attribute_wide <- as.matrix(attribute_wide)
  
  jsd_storage <- vector(mode = "numeric",length = nrow(trinuc_wide))
  
  for(jsd_ind in seq_along(jsd_storage)){
    
    jsd_storage[jsd_ind] <- suppressMessages(JSD(rbind(trinuc_wide[jsd_ind,],attribute_wide[jsd_ind,])))
    
  }
  
  return(jsd_storage)
}

# jsd_per_tumor(this_tumor = jsd_storage_list$`MDA-1006-T`)


jsd_results <- parallel::mclapply(X = jsd_storage_list,FUN = jsd_per_tumor,mc.cores = n_cores)


saveRDS(object = jsd_results, file = paste0("bootstrap_analysis/JSD_results_",tumor_name,".rds"))

purrr::map_dfr(jsd_results,median) %>%
  pivot_longer(everything(),names_to = "tumor_name",values_to = "median_JSD") %>%
  mutate(tumor_type = tumor_name) -> 
  median_jsd


saveRDS(object = median_jsd, file = paste0("bootstrap_analysis/JSD_results_median_",tumor_name,".rds"))

# jsd_results$`TCGA-O2-A5IB-01A-11D-A27K-08`
# 
# 
# jsd_results_mat <- do.call(rbind,jsd_results)

# jsd_results_mat[1:10,] %>%
#   as_tibble() %>%
#   mutate(tumor_name = rownames(jsd_results_mat)[1:10]) %>%
#   pivot_longer(cols = -tumor_name) %>%
#   ggplot(aes(x=tumor_name, y=value)) +
#   geom_jitter(height = 0,width = 0.3,alpha=0.5) +
#   theme_bw() +
#   labs(y= "JSD")


