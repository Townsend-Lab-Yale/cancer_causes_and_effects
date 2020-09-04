#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# set up data ----- 

inputs <- matrix(nrow=length(args)/2,ncol=1,data=args[seq(2,length(args),2)])
rownames(inputs) <- args[seq(1,length(args),2)]
inputs[,1] <- trimws(inputs[,1])


library(cancereffectsizeR)
analysis <- cancereffectsizeR::CESAnalysis(genome = "hg19")


# load in MAFs ---- 

if(!("-ready_MAF" %in% rownames(inputs))){
  MAF_for_analysis <- read.delim(file = inputs["-NCI_MAF",], header = T, stringsAsFactors = F,skip = inputs["-NCI_skip",])
  MAF_for_analysis$Tumor_Sample_Barcode = cancereffectsizeR::consolidate_tcga_tumors_by_patient(MAF_for_analysis$Tumor_Sample_Barcode)
  
  analysis = load_maf(analysis, maf = MAF_for_analysis, 
                      chain_file = "/ysm-gpfs/pi/townsend/general_genome_info/hg38ToHg19.over.chain")
  
  if("-Local_MAF" %in% rownames(inputs)){
    local_maf <- read.delim(file = inputs["-Local_MAF",],header = T,stringsAsFactors = F)
    
    local_maf <- dplyr::distinct(local_maf)
    if("Patient_ID" %in% colnames(local_maf)){
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
  MAF_for_analysis$Tumor_Sample_Barcode = cancereffectsizeR::consolidate_tcga_tumors_by_patient(MAF_for_analysis$Tumor_Sample_Barcode)
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

load("../mutational_signature_repertoire_clean.RData")

signatures_in_tumor_type_vec <- mutational_signature_repertoire_clean[,tumor_type_trinuc]
signatures_in_tumor_type <- rownames(mutational_signature_repertoire_clean)[which(signatures_in_tumor_type_vec == 1)]

data(signatures_names_matrix,package = "cancereffectsizeR")

signatures_to_remove <- signatures_names_matrix[!signatures_names_matrix[,1] %in% signatures_in_tumor_type,1]

# extract that artifact signatures

artifact_sigs <- signatures_names_matrix[grep(pattern = "artifact",x = signatures_names_matrix[,2]),1]

signatures_to_remove <- signatures_to_remove[!signatures_to_remove %in% artifact_sigs]


# Calculate trinucleotide mutation weightings using deconstructSigs ----
analysis <- cancereffectsizeR::trinuc_mutation_rates(analysis, cores = as.numeric(inputs["-cores",]),
                                                     signatures_to_remove = signatures_to_remove, use_dS_exome2genome=TRUE)

# Calculate gene-level mutation rates using dNdScv

analysis <- cancereffectsizeR::gene_mutation_rates(analysis, 
                                                   covariate_file = if("-cov" %in% rownames(inputs)) inputs["-cov",] else NULL)

# Assign genes to MAF, keeping assignments consistent with dndscv when possible
analysis <- cancereffectsizeR::annotate_variants(analysis)


save(analysis, file="MAF_for_analysis_after_mut_rate.RData")

analysis <-  ces_snv(cesa = analysis, 
                     cores=as.numeric(inputs["-cores",])) 

save(analysis,file = paste(inputs["-tumor_name",],"_","_selection_output.RData",sep=""))

results <- as.data.frame(analysis@selection_results[tumors_with_variant > 0])


# save(selection_output,file = paste(inputs["-tumor_name",],"_","selection_output.RData",sep=""))

save(results,file = paste(inputs["-tumor_name",],"_","_selection_output_df.RData",sep=""))



source("../population_scaled_effect_per_tumor.R")

population_scaled_selection <- population_scaled_effect_per_tumor(ces_output = analysis, min_variant_freq = 2)

source("../nrsi_postprocess.R")

population_scaled_selection <- nrsi_remove_double_annotation(analysis = analysis, nrsi_data = population_scaled_selection)

save(population_scaled_selection, file=paste(inputs["-tumor_name",],"_","population_scaled_selection.RData",sep=""))

source("../divergence_calculation.R")

tumor_JSD_df <- divergence_calculation(trinuc_weights = cancereffectsizeR::get_signature_weights(cesa = analysis,include_tumors_without_data = T), 
                                       NRSI_per_tumor = population_scaled_selection)

tumor_JSD_df$tumor_type <- inputs["-tumor_name",]

save(tumor_JSD_df, file = paste(inputs["-tumor_name",],"_","tumor_JSD_df.RData",sep=""))


q("no")
