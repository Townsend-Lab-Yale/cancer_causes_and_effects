#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# latest install with bootstrap method in mutation rate calc
# remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR")
library(cancereffectsizeR)
library(MutationalPatterns)
library(data.table)
library(tidyverse)
pbapply::pboptions(type = 'none')


source("../parameters_for_tumor_types.R")

parameters_for_run <- parameters_for_run[[ args[1] ]]

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

# thanks https://stackoverflow.com/a/53095338 
job_distribution <- split(1:n_boot, sort(rep_len(1:n_split, length(1:n_boot))))


# set.seed(1234)

# message(getwd())


# 
# 
# # preload tcga maf file 
# tcga_maf <- preload_maf(
#   maf = parameters_for_run %>%
#     dplyr::filter(input_parameter == "tcga_maf_location") %>%
#     pull(input_value),
#   refset = "ces.refset.hg19",chain_file = "/ysm-gpfs/pi/townsend/general_genome_info/hg38ToHg19.over.chain")
# tcga_maf <- tcga_maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]
# 
# 
# 
# # Create cancereffectsizeR analysis and load data
# cesa <- CESAnalysis(refset = "ces.refset.hg19")
# cesa <- load_maf(cesa = cesa, maf = tcga_maf)
# 
# # load in the YG maf if applicable
# if(
#   !is.na(parameters_for_run %>%
#          filter(input_parameter == "yg_maf_location") %>%
#          pull(input_value) ) 
# ){
#   
#   yg_maf <- preload_maf(
#     maf = parameters_for_run %>%
#       filter(input_parameter == "yg_maf_location") %>%
#       pull(input_value),
#     refset = "ces.refset.hg19",
#     chr_col = "Chrom")
#   yg_maf <- yg_maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]
#   
#   cesa <- load_maf(cesa = cesa, maf = yg_maf)  
# }
# 
# 
# 
# Infer trinucleotide-context-specific relative rates of SNV mutation from
# a mutational signature analysis (leaving out signatures not found in BRCA)
signature_exclusions <- suggest_cosmic_signatures_to_remove(cancer_type = parameters_for_run %>%
                                                              filter(input_parameter == "exclusion_cancer_type") %>%
                                                              pull(input_value),
                                                            treatment_naive = TRUE)
# 
# cesa <- gene_mutation_rates(cesa,covariates = parameters_for_run %>%
#                               filter(input_parameter == "gene_covariates") %>%
#                               pull(input_value))


cesa <- cancereffectsizeR::load_cesa(
  file = paste0("bootstrap_analysis/",args[1],"cesa_before_bootstrap.rds")
  )


n_boot <- n_boot %>%
  as.numeric()
n_cores <- n_cores %>%
  as.numeric()

trinuc_storage <- vector(mode = "list",length = n_boot)
attributable_storage <- vector(mode = "list",length = n_boot)

trinuc_storage_cesa <- vector(mode = "list",length = n_boot)

source("https://raw.githubusercontent.com/Townsend-Lab-Yale/cancer_causes_and_effects/master/R/population_scaled_effect_per_tumor.R")




these_jobs <- job_distribution[[args[2]]]



for(bootstrap_ind in these_jobs){
  
  cesa <- trinuc_mutation_rates(
    cesa = cesa, 
    signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
    signature_exclusions = signature_exclusions,
    bootstrap_mutations = TRUE,
    cores = n_cores
  )
  
  trinuc_storage[[bootstrap_ind]] <- cesa$mutational_signatures$biological_weights
  
  # trinuc_storage_cesa[[bootstrap_ind]] <- cancereffectsizeR::get_signature_weights(cesa = cesa)
  
  this_run_name <- paste0('boot_', bootstrap_ind)
  
  cesa <- ces_variant(cesa, 
                      run_name = this_run_name,
                      cores = n_cores)
  
  attributable_storage[[bootstrap_ind]] <- population_scaled_effect_per_tumor(ces_output = cesa,
                                                                              run_name = this_run_name,cores = n_cores)
  
  cesa <- clear_trinuc_rates_and_signatures(cesa)
  
  message(paste0("Bootstrap sample:", bootstrap_ind))
  
  
  # save temp files in case it breaks
  # 
  # if(bootstrap_ind %% 50 == 0){
  #   
  #   cancereffectsizeR::save_cesa(cesa = cesa,file = 
  #                                  paste0(
  #                                    "bootstrap_analysis/temp_files/",
  #                                    parameters_for_run %>%
  #                                      filter(input_parameter == "tumor_name") %>%
  #                                      pull(input_value),"_temp_cesa_",
  #                                    bootstrap_ind,
  #                                    ".rds")
  #   )
  #   saveRDS(object = trinuc_storage,
  #           file = paste0(
  #             "bootstrap_analysis/temp_files/",
  #             parameters_for_run %>%
  #               filter(input_parameter == "tumor_name") %>%
  #               pull(input_value),"_temp_trinuc_storage_",
  #             bootstrap_ind,
  #             ".rds")
  #   )
  #   saveRDS(object = attributable_storage, 
  #           file = paste0(
  #             "bootstrap_analysis/temp_files/",
  #             parameters_for_run %>%
  #               filter(input_parameter == "tumor_name") %>%
  #               pull(input_value),"_temp_attribution_storage_",
  #             bootstrap_ind,
  #             ".rds"))
  #   
  # }
  
}


if(!dir.exists("bootstrap_analysis")){
  dir.create("bootstrap_analysis")
  dir.create("bootstrap_analysis/bootstrap_results")
}

cancereffectsizeR::save_cesa(cesa = cesa,file = 
                               paste0(
                                 "bootstrap_analysis/bootstrap_results/",
                                 parameters_for_run %>%
                                   filter(input_parameter == "tumor_name") %>%
                                   pull(input_value),"_cesa_",
                                 args[2],
                                 ".rds")
)
saveRDS(object = trinuc_storage,
        file = paste0(
          "bootstrap_analysis/bootstrap_results/",
          parameters_for_run %>%
            filter(input_parameter == "tumor_name") %>%
            pull(input_value),"_trinuc_storage_",
          args[2],
          ".rds")
)
saveRDS(object = attributable_storage, 
        file = paste0(
          "bootstrap_analysis/bootstrap_results/",
          parameters_for_run %>%
            filter(input_parameter == "tumor_name") %>%
            pull(input_value),"_attribution_storage_",
          args[2],
          ".rds"))

# saveRDS(object = trinuc_storage_cesa, 
#         file = paste0(
#           "bootstrap_analysis/bootstrap_results/",
#           parameters_for_run %>%
#             filter(input_parameter == "tumor_name") %>%
#             pull(input_value),"_trinuc_storage_cesa_",
#           args[2],
#           ".rds"))


# 
# cancereffectsizeR::save_cesa(cesa = cesa,file = "~/Documents/Research/cancer_causes_and_effects/dev/bootstrap_tests/cesa_boot_test_LUSC.rds")
# saveRDS(object = trinuc_storage,file = "~/Documents/Research/cancer_causes_and_effects/dev/bootstrap_tests/trinuc_storage_boot_test_LUSC.rds")
# saveRDS(object = attributable_storage, file = "~/Documents/Research/cancer_causes_and_effects/dev/bootstrap_tests/attributable_storage_LUSC.rds")
