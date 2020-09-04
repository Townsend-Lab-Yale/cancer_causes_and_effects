library(tidyverse)

# Selecting one annotation for each variant ---- 
nrsi_remove_double_annotation <- function(analysis, nrsi_data){
  new_nrsi <- nrsi_data
  message("Selecting one annotation for each variant...")
  for(tumor_ind in seq_along(unique(new_nrsi$Unique_Patient_Identifier))){

    # if(tumor_ind==60){break}
    
    # subset the data of interest 
    this_tumor <- unique(new_nrsi$Unique_Patient_Identifier)[tumor_ind]
    this_maf <- analysis@maf %>%
      filter(Unique_Patient_Identifier == this_tumor)
    this_nrsi <- new_nrsi %>%
      filter(Unique_Patient_Identifier == this_tumor)
    
    # see if there are any double annotations that are recurrent 
    # and in the nrsi output
    if(any(sapply(this_maf$assoc_aa_mut,length) > 1)){
      double_attribute <- 
        unlist(this_maf[which(sapply(this_maf$assoc_aa_mut,length) > 1),
                        assoc_aa_mut])
      if(any(double_attribute %in% this_nrsi$variant)){
        # double attributes are in the nrsi data. Need to remove one
        # choose to remove the attribute with the lowest fitness (nrsi) 
        double_attribute <- this_maf[which(sapply(this_maf$assoc_aa_mut,length) > 1),
                                     assoc_aa_mut]
        for(attributes in seq_along(double_attribute)){
          
          these_attributes <- double_attribute[[attributes]]
          
          
          if(length(which(these_attributes %in% this_nrsi$variant))>1){
            total_nrsi <- this_nrsi %>% 
              filter(variant %in% these_attributes) %>%
              dplyr::select(ends_with("nrsi")) %>% 
              rowSums()
            
            variant_to_keep <- these_attributes[which.max(total_nrsi)]
            variants_to_drop <- these_attributes[-which.max(total_nrsi)]
            
            new_nrsi <- new_nrsi[-which(
              new_nrsi$variant %in% variants_to_drop & 
                new_nrsi$Unique_Patient_Identifier == this_tumor),]
            
          }
          
        }
        
      }
      
    }
  }
  return(new_nrsi)
}
