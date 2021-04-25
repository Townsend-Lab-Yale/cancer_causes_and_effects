#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# set up data ----- 

inputs <- matrix(nrow=length(args)/2,ncol=1,data=args[seq(2,length(args),2)])
rownames(inputs) <- args[seq(1,length(args),2)]
inputs[,1] <- trimws(inputs[,1])

library(data.table)
library(cancereffectsizeR)
analysis <- cancereffectsizeR::CESAnalysis(refset = "ces.refset.hg19")


# load in MAFs ---- 

if(!("-ready_MAF" %in% rownames(inputs))){
  MAF_for_analysis <- read.delim(file = inputs["-NCI_MAF",], header = T, stringsAsFactors = F,skip = inputs["-NCI_skip",])
  MAF_for_analysis$Tumor_Sample_Barcode = cancereffectsizeR:::consolidate_tcga_tumors_by_patient(MAF_for_analysis$Tumor_Sample_Barcode)
  
  analysis = load_maf(analysis, maf = MAF_for_analysis, 
                      chain_file = "/ysm-gpfs/pi/townsend/general_genome_info/hg38ToHg19.over.chain")
  
  if("-Local_MAF" %in% rownames(inputs)){
    local_maf <- read.delim(file = inputs["-Local_MAF",],header = T,stringsAsFactors = F)
    
    local_maf <- dplyr::distinct(local_maf)
    if("Patient_ID" %in% colnames(local_maf) &
       !("Tumor_Sample_Barcode" %in% colnames(local_maf))){
      colnames(local_maf)[which(colnames(local_maf) == "Patient_ID") ] <- "Tumor_Sample_Barcode"
    }
    
    
    if("Chrom" %in% colnames(local_maf)){
      colnames(local_maf)[which(colnames(local_maf) == "Chrom") ] <- "Chromosome"
    }
    
    # MAF_for_analysis <- cancereffectsizeR::merging_TCGA_and_local_MAFdata_function(NCI_data = MAF_for_analysis,Local_data = local_maf,check_for_same_tumor = T)
    analysis <- load_maf(analysis, maf = local_maf)
  }
  
}else{
  MAF_for_analysis <- get(load(inputs["-ready_MAF",]))
  
  MAF_for_analysis <- dplyr::distinct(MAF_for_analysis[,c("Chromosome","Start_Position","Tumor_Sample_Barcode","Tumor_Seq_Allele2","Tumor_Seq_Allele1","Reference_Allele")])
  MAF_for_analysis$Tumor_Sample_Barcode = cancereffectsizeR:::consolidate_tcga_tumors_by_patient(MAF_for_analysis$Tumor_Sample_Barcode)
  analysis <- load_maf(analysis, maf = MAF_for_analysis)
  
  
}


# remove any YG tumors that had therapy prior to sequencing ----
load("../YG_tumors_with_therapy.RData")
if(any(analysis@maf$Unique_Patient_Identifier %in% tumors_with_therapy)){
  removed_tumors <- analysis@maf$Unique_Patient_Identifier[which(analysis@maf$Unique_Patient_Identifier %in% tumors_with_therapy)]
  save(removed_tumors,file = paste(inputs["-tumor_name",],"_","removed_tumors.RData",sep=""))
  analysis@maf <- analysis@maf[-which(analysis@maf$Unique_Patient_Identifier %in% tumors_with_therapy),]
}


MAF_for_analysis <- analysis@maf

save(MAF_for_analysis, file = paste(inputs["-tumor_name",],"_","MAF.RData",sep=""))

tumor_type_trinuc <- as.character(unlist(strsplit(inputs["-tumor_type_trinuc_sig",], ",")))


# ESCA is actually made up of ESCA and ESCC, which 
# have different possible trinucleotide signatures... 
# breaking them down here and putting them back together
if(inputs["-tumor_name",] == "ESCA"){
  
  MAF_for_analysis$tumor_esca_type <- ""
  
  esca_tumor_list <- read.csv(file = "our_esca_tumors.csv")
  esca_tumor_list$CASE_ID <- substr(esca_tumor_list$CASE_ID,start = 1,stop = 12)
  
  
  for(tumor_ind in seq_along(unique(MAF_for_analysis$Unique_Patient_Identifier))){
    MAF_for_analysis[MAF_for_analysis$Unique_Patient_Identifier == unique(MAF_for_analysis$Unique_Patient_Identifier)[tumor_ind],"tumor_esca_type"] <- 
      esca_tumor_list[esca_tumor_list$CASE_ID == unique(MAF_for_analysis$Unique_Patient_Identifier)[tumor_ind],"ONCOTREE_CODE"]
  }
  
  
  analysis <- cancereffectsizeR::CESAnalysis(refset = "ces.refset.hg19",sample_groups = c("ESCA","ESCC"))
  analysis <- cancereffectsizeR::load_maf(cesa = analysis, 
                                          maf = MAF_for_analysis, 
                                          group_col = "tumor_esca_type")
  
  # trinucs on ESCA
  analysis <- cancereffectsizeR::trinuc_mutation_rates(cesa = analysis, 
                                                       signature_set = "COSMIC_v3.1", 
                                                       cores = as.numeric(inputs["-cores",]),
                                                       signatures_to_remove = 
                                                         c(as.character(cancereffectsizeR::suggest_cosmic_signatures_to_remove(
                                                           cancer_type = "Eso-AdenoCA",treatment_naive = T
                                                         )),"SBS89"), 
                                                       use_dS_exome2genome=TRUE,sample_group = "ESCA")
  
  # trinucs on ESCC
  # include signature 16 in accordance with https://www.annalsofoncology.org/article/S0923-7534(19)45461-6/fulltext
  
  escc_to_remove <- c(as.character(cancereffectsizeR::suggest_cosmic_signatures_to_remove(
    cancer_type = "Eso-SCC",treatment_naive = T
  )),"SBS89")
  
  escc_to_remove <- escc_to_remove[-which(escc_to_remove == "SBS16")]
  
  analysis <- cancereffectsizeR::trinuc_mutation_rates(cesa = analysis, 
                                                       signature_set = "COSMIC_v3.1", 
                                                       cores = as.numeric(inputs["-cores",]),
                                                       signatures_to_remove = escc_to_remove, 
                                                       use_dS_exome2genome=TRUE,sample_group = "ESCC")
  
  
}else{
  
  # signatures_names_matrix_new <- cancereffectsizeR::get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.1")
  # signatures_names_matrix_new <- signatures_names_matrix_new$meta[,c("Signature","Etiology")]
  to_remove = c("SBS89",as.character(suggest_cosmic_signatures_to_remove(cancer_type = tumor_type_trinuc, treatment_naive = T)))
  
  if(tumor_type_trinuc %in% c("HNSC","BLCA")){
    if("SBS88" %in% to_remove){
      to_remove <- to_remove[-which(to_remove == "SBS88")]
    }
  }
  
  # Calculate trinucleotide mutation weightings using deconstructSigs ----
  analysis <- cancereffectsizeR::trinuc_mutation_rates(analysis,signature_set = "COSMIC_v3.1", cores = as.numeric(inputs["-cores",]),
                                                       signatures_to_remove = to_remove, use_dS_exome2genome=TRUE)
  
  
  
}



# Calculate gene-level mutation rates using dNdScv

analysis <- cancereffectsizeR::gene_mutation_rates(analysis, 
                                                   covariates = if("-cov" %in% rownames(inputs)) inputs["-cov",] else NULL)



save(analysis, file="MAF_for_analysis_after_mut_rate.RData")

analysis <-  ces_variant(cesa = analysis, 
                         cores=as.numeric(inputs["-cores",])) 

save(analysis,file = paste(inputs["-tumor_name",],"_","_selection_output.RData",sep=""))

# results <- as.data.frame(analysis@selection_results[tumors_with_variant > 0])


# save(selection_output,file = paste(inputs["-tumor_name",],"_","selection_output.RData",sep=""))

# save(results,file = paste(inputs["-tumor_name",],"_","_selection_output_df.RData",sep=""))



source("../population_scaled_effect_per_tumor.R")

population_scaled_selection <- population_scaled_effect_per_tumor(ces_output = analysis, min_variant_freq = 2)

source("../nrsi_postprocess.R")

population_scaled_selection <- nrsi_remove_double_annotation(analysis = analysis, nrsi_data = population_scaled_selection)

save(population_scaled_selection, file=paste(inputs["-tumor_name",],"_","population_scaled_selection.RData",sep=""))

source("../divergence_calculation.R")

tumor_JSD_df <- divergence_calculation(trinuc_weights = cancereffectsizeR::get_signature_weights(cesa = analysis), 
                                       NRSI_per_tumor = population_scaled_selection)

tumor_JSD_df$tumor_type <- inputs["-tumor_name",]

save(tumor_JSD_df, file = paste(inputs["-tumor_name",],"_","tumor_JSD_df.RData",sep=""))


q("no")
