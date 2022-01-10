
library(cancereffectsizeR)
library(data.table)
library(tidyverse)

files <- dir(recursive = T,full.names = T,pattern = "cesa_before")

variant_number_extractor <- function(this_file){
  this_cesa <- readRDS(this_file)
  
  variants <- this_cesa$maf %>%
    count(Unique_Patient_Identifier)
  
  tumor_name <- strsplit(this_file,"cesa")[[1]]
  tumor_name <- tumor_name[1]
  tumor_name <- strsplit(tumor_name,"/")[[1]]
  tumor_name <- tumor_name[length(tumor_name)]
  
  variants <- variants %>% mutate(tumor_name = tumor_name)
  
  return(variants)
  
  
}

variant_numbers <- purrr::map_dfr(files, variant_number_extractor)

saveRDS(object = variant_numbers, file = "variant_numbers.rds") 