
#' Some ESCA TCGA tumors are actually ESCC (squamous cell), which have different
#' detected trinucleotide signatures than ESCA. This script uses
#' the TCGAretriever package to get the clinical data for ESCA tumors, and 
#' parse out which of our tumors belongs to ESCA or ESCC. Then, in the 
#' CESAnalysis code, we will apply tumor-specific trinucleotide calculations
#' depending on whether the `ESCA` tumor is ESCA or ESCC. 


library(TCGAretriever)


tcga_esca_clin <- TCGAretriever::get_clinical_data(case_id = "esca_tcga_all") # find ESCA clinical data
tcga_esca_clin_code <- tcga_esca_clin[,c("CASE_ID","ONCOTREE_CODE")] # extract just the tumor type code

# our ESCA MAF file from TCGA
tcga_esca <- read.delim(file = "dev/ESCA/TCGA.ESCA.mutect.7f8e1e7c-621c-4dfd-8fad-af07c739dbfc.DR-10.0.somatic.maf",header = T,skip = 5) 

# extract our tumor names 
esca_tumors <- unique(tcga_esca$Tumor_Sample_Barcode)
esca_tumors <- substr(x = esca_tumors,start = 1,stop = 15)


rownames(tcga_esca_clin_code) <- tcga_esca_clin_code$CASE_ID
esca_tumors %in% tcga_esca_clin_code$CASE_ID # clinical data has all our cases! 

table(tcga_esca_clin_code[esca_tumors,"ONCOTREE_CODE"])

# extract just our tumors and save for use on the cluster
our_esca_tumors <- tcga_esca_clin_code[esca_tumors,]
# save(our_esca_tumors,file = "dev/our_esca_tumors.RData")

write.table(x = our_esca_tumors,file = "dev/our_esca_tumors.csv",quote = F,row.names = F,sep = ",")

# esca_read_test <- read.csv(file = "dev/our_esca_tumors.csv")


tcga_esca_clin_code

tcga_esca$short_barcode <- substr(x = tcga_esca$Tumor_Sample_Barcode,start = 1,stop = 15)

tcga_esca_esca <- tcga_esca[tcga_esca$short_barcode %in% tcga_esca_clin_code$CASE_ID[tcga_esca_clin_code$ONCOTREE_CODE=="ESCA"],]
tcga_esca_escc <- tcga_esca[tcga_esca$short_barcode %in% tcga_esca_clin_code$CASE_ID[tcga_esca_clin_code$ONCOTREE_CODE=="ESCC"],]

write_delim(x = tcga_esca_esca, file = "dev/ESCA/TCGA_ESCA_just_esca.tsv")
write_delim(x = tcga_esca_escc, file = "dev/ESCA/TCGA_ESCA_just_escc.tsv")



tcga_esca <- preload_maf(maf = tcga_esca_esca,refset = "ces.refset.hg19",chain_file = "dev/hg38ToHg19.over.chain")
tcga_escc <- preload_maf(maf = tcga_esca_escc,refset = "ces.refset.hg19",chain_file = "dev/hg38ToHg19.over.chain")


saveRDS(object = tcga_esca, file = "dev/combine_test/ESCA/tcga_esca.rds")
saveRDS(object = tcga_escc, file = "dev/combine_test/ESCA/tcga_escc.rds")


save(object = tcga_esca, file = "dev/combine_test/ESCA/tcga_esca.RData")
save(object = tcga_escc, file = "dev/combine_test/ESCA/tcga_escc.RData")


