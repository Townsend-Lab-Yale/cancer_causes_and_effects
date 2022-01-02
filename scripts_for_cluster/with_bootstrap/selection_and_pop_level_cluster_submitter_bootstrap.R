#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# 
# inputs <- matrix(nrow=length(args)/2,ncol=1,data=args[seq(2,length(args),2)])
# rownames(inputs) <- args[seq(1,length(args),2)]
# inputs[,1] <- trimws(inputs[,1])
# 
library(tidyverse)
library(cancereffectsizeR)


source("../parameters_for_tumor_types.R")

parameters_for_run <- parameters_for_run[[args]]

n_cores <- parameters_for_run %>% 
  filter(input_parameter == "n_cores") %>% 
  pull(input_value)

n_boot <- parameters_for_run %>% 
  filter(input_parameter == "n_boot") %>% 
  pull(input_value)

n_split <- parameters_for_run %>% 
  filter(input_parameter == "split") %>% 
  pull(input_value)


tumor_name <- parameters_for_run %>% 
  filter(input_parameter == "tumor_name") %>% 
  pull(input_value)


partition <- parameters_for_run %>% 
  filter(input_parameter == "partition") %>% 
  pull(input_value)


time <- parameters_for_run %>% 
  filter(input_parameter == "time") %>% 
  pull(input_value)

memory <- parameters_for_run %>% 
  filter(input_parameter == "memory") %>% 
  pull(input_value)

# thanks https://stackoverflow.com/a/53095338 
job_distribution <- split(1:n_boot, sort(rep_len(1:n_split, length(1:n_boot))))



# Create cancereffectsizeR analysis and load data ---- 
cesa <- CESAnalysis(refset = "ces.refset.hg19")



if(
  !is.na(parameters_for_run %>%
         dplyr::filter(input_parameter == "tcga_maf_location") %>%
         pull(input_value))
){
  # preload tcga maf file 
  tcga_maf <- preload_maf(
    maf = parameters_for_run %>%
      dplyr::filter(input_parameter == "tcga_maf_location") %>%
      pull(input_value),
    refset = "ces.refset.hg19",chain_file = "/ysm-gpfs/pi/townsend/general_genome_info/hg38ToHg19.over.chain")
  tcga_maf <- tcga_maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]
  
  
  cesa <- load_maf(cesa = cesa, maf = tcga_maf)
  
}



# load in the YG maf if applicable
if(
  !is.na(parameters_for_run %>%
         filter(input_parameter == "yg_maf_location") %>%
         pull(input_value) ) 
){
  
  yg_maf <- preload_maf(
    maf = parameters_for_run %>%
      filter(input_parameter == "yg_maf_location") %>%
      pull(input_value),
    refset = "ces.refset.hg19",
    chr_col = "Chrom")
  yg_maf <- yg_maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]
  
  cesa <- load_maf(cesa = cesa, maf = yg_maf)  
}



if(
  !is.na(parameters_for_run %>%
         dplyr::filter(input_parameter == "ready_maf") %>%
         pull(input_value))
){
  
  MAF_for_analysis <- get(
    load(
      parameters_for_run %>%
        dplyr::filter(input_parameter == "ready_maf") %>%
        pull(input_value)
    )
  )
  
  if(!tumor_name %in% c("ESCC","ESCA")){
    MAF_for_analysis <- dplyr::distinct(MAF_for_analysis[,c("Chromosome","Start_Position","Tumor_Sample_Barcode","Tumor_Seq_Allele2","Tumor_Seq_Allele1","Reference_Allele")])
  }
  cesa <- load_maf(cesa = cesa, maf = MAF_for_analysis)
  
  
}





# remove any YG tumors that had therapy prior to sequencing ----
load("../YG_tumors_with_therapy.RData")
if(any(cesa@maf$Unique_Patient_Identifier %in% tumors_with_therapy)){
  removed_tumors <- cesa@maf$Unique_Patient_Identifier[which(cesa@maf$Unique_Patient_Identifier %in% tumors_with_therapy)]
  save(removed_tumors,file = paste(parameters_for_run %>%
                                     filter(input_parameter == "tumor_name") %>%
                                     pull(input_value),"_","removed_tumors.RData",sep=""))
  cesa@maf <- cesa@maf[-which(cesa@maf$Unique_Patient_Identifier %in% tumors_with_therapy),]
}



# Infer trinucleotide-context-specific relative rates of SNV mutation from
# a mutational signature analysis (leaving out signatures not found in BRCA)
signature_exclusions <- suggest_cosmic_signatures_to_remove(cancer_type = parameters_for_run %>%
                                                              filter(input_parameter == "exclusion_cancer_type") %>%
                                                              pull(input_value), 
                                                            treatment_naive = TRUE)

cesa <- gene_mutation_rates(cesa,covariates = parameters_for_run %>%
                              filter(input_parameter == "gene_covariates") %>%
                              pull(input_value))

if(!dir.exists("bootstrap_analysis")){
  dir.create("bootstrap_analysis")
  dir.create("bootstrap_analysis/bootstrap_results")
}

cancereffectsizeR::save_cesa(cesa = cesa, 
                             file = paste0("bootstrap_analysis/",args,"cesa_before_bootstrap.rds")
)




# Create batch job and submit to the cluster for each bootstrap ----

for(job_ind in seq_along(job_distribution)){
  job_name <- paste0(tumor_name,"_",job_ind)  
  
  file.create(paste("selection_job_",args,"_",job_ind,".sh",sep=""))
  sink(file = paste("selection_job_",args,"_",job_ind,".sh",sep=""))
  cat("#!/bin/bash\n")
  cat(paste("#SBATCH -n 1\n#SBATCH -c ",n_cores,"\n#SBATCH -J ", job_name ,"_selection\n#SBATCH --partition=",partition,"\n#SBATCH -t ", time ,"\n#SBATCH --mem-per-cpu=",memory,"\n#SBATCH --mail-user=vincent.cannataro@yale.edu\n#SBATCH --mail-type=FAIL\n",sep=""),append=T)
  cat(paste("module load GCC/10.2.0\n"),append = T)
  cat(paste("~/../../pi/townsend/vlc24/R_source/R-4.1.2/bin/Rscript ../selection_and_pop_level_cluster_bootstrap.R",args,job_ind,sep=" "),append = T)
  sink()
  system(paste("sbatch selection_job_",args,"_",job_ind,".sh",sep=""))
  
}


# run the no bootstrap case ----
source("https://raw.githubusercontent.com/Townsend-Lab-Yale/cancer_causes_and_effects/master/R/population_scaled_effect_per_tumor.R")



cesa <- trinuc_mutation_rates(
  cesa = cesa, 
  signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
  signature_exclusions = signature_exclusions,
  bootstrap_mutations = FALSE,
  cores = n_cores
)


trinuc_storage <- cesa$mutational_signatures$biological_weights

trinuc_storage_cesa <- cancereffectsizeR::get_signature_weights(cesa = cesa)


cesa <- ces_variant(cesa, 
                    run_name = "no_bootstrap",
                    cores = n_cores)

attributable_storage <- population_scaled_effect_per_tumor(ces_output = cesa,
                                                           run_name = "no_bootstrap",
                                                           cores = n_cores)




if(!dir.exists("no_bootstrap_analysis")){
  dir.create("no_bootstrap_analysis")
  # dir.create("bootstrap_analysis/bootstrap_results")
}

cancereffectsizeR::save_cesa(cesa = cesa,file = 
                               paste0(
                                 "no_bootstrap_analysis/",
                                 parameters_for_run %>%
                                   filter(input_parameter == "tumor_name") %>%
                                   pull(input_value),"_cesa",
                                 ".rds")
)
saveRDS(object = trinuc_storage,
        file = paste0(
          "no_bootstrap_analysis/",
          parameters_for_run %>%
            filter(input_parameter == "tumor_name") %>%
            pull(input_value),"_trinuc_storage",
          ".rds")
)
saveRDS(object = attributable_storage, 
        file = paste0(
          "no_bootstrap_analysis/",
          parameters_for_run %>%
            filter(input_parameter == "tumor_name") %>%
            pull(input_value),"_attribution_storage",
          ".rds"))

saveRDS(object = trinuc_storage_cesa, 
        file = paste0(
          "no_bootstrap_analysis/",
          parameters_for_run %>%
            filter(input_parameter == "tumor_name") %>%
            pull(input_value),"_trinuc_storage_cesa",
          ".rds"))




