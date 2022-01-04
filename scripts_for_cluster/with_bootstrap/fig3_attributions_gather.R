

library(tidyverse)
library(data.table)

# the following file is from this manuscript:
# Bailey, Matthew H., Collin Tokheim, Eduard Porta-Pardo, Sohini Sengupta, Denis Bertrand, Amila Weerasinghe, Antonio Colaprico, et al. 2018. “Comprehensive Characterization of Cancer Driver Genes and Mutations.” Cell 174 (4): 1034–35.

# The file is located here: https://www.cell.com/cms/10.1016/j.cell.2018.02.060/attachment/b0abfafa-a1ff-4f7d-8842-b49ba8d32e08/mmc1.xlsx

# I will not include the file in our repo as I am unsure if we have permission to
# share/host. But, it is open access!
bailey_driver_list <- readxl::read_excel(path = "Bailey_etal_Cell_1-s2.0-S009286741830237X-mmc1.xlsx",sheet = "Table S1",skip = 3)


tumors_for_attribution_fig3 <- c("LUAD","LUSC","SKCMP","LIHC","BLCA","CESC","HNSC_HPVneg","HNSC_HPVpos")

file_names <- paste0("combined_attribution_",tumors_for_attribution_fig3)

files_names <- dir(recursive = T,pattern = paste0(file_names,collapse = "|"),full.names = T)

# files_names <- files_names[2:3] # for dev


attribution_storage <- vector(mode = "list",length = length(tumors_for_attribution_fig3))
names(attribution_storage) <- tumors_for_attribution_fig3

trinuc_storage <- vector(mode = "list",length = length(tumors_for_attribution_fig3))
names(trinuc_storage) <- tumors_for_attribution_fig3


for( file_ind in seq_along(files_names)){
  
  these_attributions <- readRDS(file = files_names[file_ind])
  these_attributions <- rbindlist(these_attributions, idcol = "bootstrap_sample",fill = T)
  
  # different tumors had different attribution columns, fill in NA above then replace with zero (there was zero attribution to that missing signature if missing)
  these_attributions <- these_attributions %>%
    replace(is.na(.),0)
  
  
  
  num_of_tumors <- length(unique(these_attributions$Unique_Patient_Identifier))
  
  attributions_all_variants <- these_attributions %>%
    select(-ends_with("flux")) %>%
    pivot_longer(cols = starts_with("SBS"),names_to = "signature",values_to = "weight") %>%
    filter(weight > 0) %>%
    group_by(bootstrap_sample,Unique_Patient_Identifier) %>%
    mutate(proportional_weight = weight / sum(weight)) %>%
    group_by(bootstrap_sample,variant,signature) %>%
    summarize(avg_weight = sum(proportional_weight)/num_of_tumors) %>% 
    ungroup()
    
  attributions_all_variants %>% 
    group_by(variant,signature) %>%
    summarize(median_weight = median(avg_weight)) %>%
    group_by(variant) %>% 
    summarize(total_med = sum(median_weight)) %>%
    arrange(desc(total_med)) %>%
    mutate(gene = stringr::word(string = variant,start = 1,sep = "_")) -> 
    top_attribution
    
  
  
  this_file <- files_names[file_ind]
  subnames <- str_split(this_file,pattern = ".rds")[[1]]
  subnames <- str_split(subnames[1],pattern = "combined_attribution_")[[1]]
  tumor_type <- subnames[length(subnames)]
  bailey_tumor_type <- tumor_type
  
  if(bailey_tumor_type == "HNSC_HPVneg"){
    bailey_tumor_type <- "HNSC"
  }
  
  if(bailey_tumor_type == "HNSC_HPVpos"){
    bailey_tumor_type <- "HNSC"
  }
  
  if(bailey_tumor_type == "SKCMP"){
    bailey_tumor_type <- "SKCM"
  }
  
  if(bailey_tumor_type == "SKCMM"){
    bailey_tumor_type <- "SKCM"
  }

  
  bailey_driver_list_subset <- bailey_driver_list %>% 
    filter(Cancer %in% c("PANCAN",bailey_tumor_type))
  
  top_attribution_variants <- top_attribution %>%
    filter(gene %in% bailey_driver_list_subset$Gene) %>%
    pull(variant)
  
  
  top_attribution_variants <- top_attribution_variants[1:ifelse(length(top_attribution_variants)>=10,10,length(top_attribution_variants))]

  
  
  
  attribution_top_variants <- attributions_all_variants %>% 
    filter(variant %in% top_attribution_variants)

  # word(attribution_top_variants$variant,1,2,"_")
  
  
  
  
  attribution_storage[[tumor_type]] <- attribution_top_variants
  
  
  # trinuc analysis for avg in among bootstraps
  
  
  trinuc_file <- stringr::str_replace(string = files_names[file_ind],
                                      pattern = "combined_attribution_",
                                      replacement = "combined_trinuc_storage_")
  
  these_trinuc <- readRDS(file = trinuc_file)
  these_trinuc <- rbindlist(these_trinuc,idcol="bootstrap_sample")
  
  these_trinuc %>% 
    select(-ends_with("snvs"),-group_avg_blended) %>%
    pivot_longer(cols = starts_with("SBS"),names_to = "signature",values_to = "weight") %>%
    group_by(signature) %>% 
    summarize(avg_weight = mean(weight)) %>% 
    arrange(desc(avg_weight)) -> 
    trinuc_avg_weight
  
  
  
  trinuc_storage[[tumor_type]] <- trinuc_avg_weight
  
  print(file_ind)
  
  
  
}


saveRDS(object = attribution_storage, file = "fig3_data_gathered.rds")
saveRDS(object = trinuc_storage, file = "fig3_data_trinuc_gathered.rds")
