

parameters_for_run <- list()

library(tibble)


parameters_for_run[["UCEC"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.UCEC.mutect.d3fa70be-520a-420e-bb6d-651aeee5cb50.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "UCEC",
  "gene_covariates",       "UCEC",
  "n_boot",                "1000",
  "n_cores",               "3",
  "split",                 "100",
  "tumor_name",            "UCEC",
  "partition",             "pi_townsend",
  "time",                  "10:00:00",
  "memory",                "20000",
  "ready_maf",             NA
)

parameters_for_run[["UCEC_bootcombine"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.UCEC.mutect.d3fa70be-520a-420e-bb6d-651aeee5cb50.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "UCEC",
  "gene_covariates",       "UCEC",
  "n_boot",                "1000",
  "n_cores",               "3",
  "split",                 "100",
  "tumor_name",            "UCEC",
  "partition",             "pi_townsend",
  "time",                  "10:00:00",
  "memory",                "30000",
  "ready_maf",             NA
)



parameters_for_run[["LUSC"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.LUSC.mutect.95258183-63ea-4c97-ae29-1bae9ed06334.DR-10.0.somatic.maf",
  "yg_maf_location",       "mutationsTN_63_Lung_Squamous_Cell_Carcinoma_Yale___MD_Anderson.maf",
  "exclusion_cancer_type", "LUSC",
  "gene_covariates",       "lung",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "LUSC",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)



parameters_for_run[["BLCA"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.BLCA.mutect.0e239d8f-47b0-4e47-9716-e9ecc87605b9.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "BLCA",
  "gene_covariates",       "bladder",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "BLCA",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)  


parameters_for_run[["BRCA_ER_neg"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "BRCA_NCI_ER_neg.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "BRCA",
  "gene_covariates",       "breast",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "BRCA_ER_neg",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)  






parameters_for_run[["BRCA_ER_pos"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "BRCA_NCI_ER_pos.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "BRCA",
  "gene_covariates",       "breast",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "BRCA_ER_pos",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)  




parameters_for_run[["CESC"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.CESC.mutect.5ffa70b1-61b4-43d1-b10a-eda412187c17.DR-10.0.somatic.maf",
  "yg_maf_location",       "mutationsTN_32_Cervical_cancer_no_cell_lines_just_SCC.maf",
  "exclusion_cancer_type", "CESC",
  "gene_covariates",       "CESC",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "CESC",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)  


parameters_for_run[["COAD"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf",
  "yg_maf_location",       "mutationsTN_41_Metastatic_colon_cancer_to_liver__lung_and_elsewhere_JUST_primary.maf",
  "exclusion_cancer_type", "COAD",
  "gene_covariates",       "colon",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "COAD",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "15000",
  "ready_maf",             NA
)  


#TODO sort out ESCA and ESCC

parameters_for_run[["ESCA"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     NA,
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "Eso-AdenoCA",
  "gene_covariates",       "ESCA",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "ESCA",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             "tcga_esca.RData"
)  


parameters_for_run[["ESCC"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     NA,
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "Eso-SCC",
  "gene_covariates",       "ESCA",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "ESCC",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             "tcga_escc.RData"
)  


parameters_for_run[["GBM"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.GBM.mutect.da904cd3-79d7-4ae3-b6c0-e7127998b3e6.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "GBM",
  "gene_covariates",       "GBM",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "GBM",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)  




parameters_for_run[["HNSC_HPVneg"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     NA,
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "HNSC",
  "gene_covariates",       "HNSC",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "HNSC_HPVneg",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             "HNSC_HPV_neg.RData"
)  






parameters_for_run[["HNSC_HPVpos"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     NA,
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "HNSC",
  "gene_covariates",       "HNSC",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "HNSC_HPVpos",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             "HNSC_HPV_pos.RData"
)  

# 
# 


parameters_for_run[["KIRC"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.KIRC.mutect.2a8f2c83-8b5e-4987-8dbf-01f7ee24dc26.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "KIRC",
  "gene_covariates",       "KIRC",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "KIRC",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)  



parameters_for_run[["LGG"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.LGG.mutect.1e0694ca-fcde-41d3-9ae3-47cfaf527f25.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "LGG",
  "gene_covariates",       "LGG",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "LGG",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)  



parameters_for_run[["LUAD"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf",
  "yg_maf_location",       "mutationsTN_62_Lung_Adenocarcinoma_Yale___MD_Anderson.maf",
  "exclusion_cancer_type", "LUAD",
  "gene_covariates",       "lung",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "LUAD",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)  

 

parameters_for_run[["LIHC"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.LIHC.mutect.a630f0a0-39b3-4aab-8181-89c1dde8d3e2.DR-10.0.somatic.maf",
  "yg_maf_location",       "mutationsTN_29_Hepatocellular_Carcinoma.maf",
  "exclusion_cancer_type", "LIHC",
  "gene_covariates",       "LIHC",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "LIHC",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "5000",
  "ready_maf",             NA
)  


parameters_for_run[["OV"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.OV.mutect.b22b85eb-2ca8-4c9f-a1cd-b77caab999bd.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "OV",
  "gene_covariates",       "OV",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "OV",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "5000",
  "ready_maf",             NA
)  


parameters_for_run[["PAAD"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf",
  "yg_maf_location",       "new_mutationsTN_26_Pancreatic_Cancer.maf",
  "exclusion_cancer_type", "PAAD",
  "gene_covariates",       "pancreas",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "PAAD",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "5000",
  "ready_maf",             NA
)  





parameters_for_run[["PRAD"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.PRAD.mutect.deca36be-bf05-441a-b2e4-394228f23fbe.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "PRAD",
  "gene_covariates",       "PRAD",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "PRAD",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "5000",
  "ready_maf",             NA
)  


parameters_for_run[["READ"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "READ",
  "gene_covariates",       "rectum",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "READ",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "5000",
  "ready_maf",             NA
)  


parameters_for_run[["SKCM_primary"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     NA,
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "SKCM",
  "gene_covariates",       "SKCM",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "SKCMP",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "15000",
  "ready_maf",             "SKCM_NCI_Yale_prim.RData"
)  


parameters_for_run[["SKCM_metastasis"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     NA,
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "SKCM",
  "gene_covariates",       "SKCM",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "SKCMM",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "15000",
  "ready_maf",             "SKCM_NCI_Yale_met.RData"
)  


parameters_for_run[["STAD"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.STAD.mutect.c06465a3-50e7-46f7-b2dd-7bd654ca206b.DR-10.0.somatic.maf",
  "yg_maf_location",       NA,
  "exclusion_cancer_type", "Stomach-AdenoCA",
  "gene_covariates",       "stomach",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "STAD",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "10000",
  "ready_maf",             NA
)  




parameters_for_run[["THCA"]] <- tibble::tribble(
  ~input_parameter,     ~input_value,
  # -----------------------|--------------
  "tcga_maf_location",     "TCGA.THCA.mutect.13999735-2e70-439f-a6d9-45d831ba1a1a.DR-10.0.somatic.maf",
  "yg_maf_location",       "mutationsTN_36_Anaplastic_thyroid_carcinoma.maf",
  "exclusion_cancer_type", "THCA",
  "gene_covariates",       "THCA",
  "n_boot",                "1000",
  "n_cores",               "5",
  "split",                 "10",
  "tumor_name",            "THCA",
  "partition",             "general",
  "time",                  "10:00:00",
  "memory",                "8000",
  "ready_maf",             NA
)  
#   