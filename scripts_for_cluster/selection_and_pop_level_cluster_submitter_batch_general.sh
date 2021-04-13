#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -J ALL_selection
#SBATCH --partition=general
#SBATCH -C haswell
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-user=vincent.cannataro@yale.edu
#SBATCH --mail-type=FAIL
cd UCEC; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R  -NCI_MAF TCGA.UCEC.mutect.d3fa70be-520a-420e-bb6d-651aeee5cb50.DR-10.0.somatic.maf -tumor_name UCEC -partition general -cov_name UCEC -NCI_skip 5 -cores 8 -mem 14000  -time 120:00:00 -tumor_type_trinuc_sig UCEC; cd ../
cd BLCA; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.BLCA.mutect.0e239d8f-47b0-4e47-9716-e9ecc87605b9.DR-10.0.somatic.maf -partition general -tumor_name BLCA -NCI_skip 5 -cov bladder -cores 12 -mem 5300 -time 24:00:00 -tumor_type_trinuc_sig BLCA; cd ../
cd BRCA_ER_neg; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF BRCA_NCI_ER_neg.maf -partition general -tumor_name BRCA_ER_neg -NCI_skip 0 -cov breast -cores 12 -mem 5300 -time 24:00:00 -tumor_type_trinuc_sig BRCA; cd ../
cd BRCA_ER_pos; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF BRCA_NCI_ER_pos.maf -partition general -tumor_name BRCA_ER_pos -NCI_skip 0 -cov breast -cores 12 -mem 5300 -time 24:00:00 -tumor_type_trinuc_sig BRCA; cd ../
cd CESC; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.CESC.mutect.5ffa70b1-61b4-43d1-b10a-eda412187c17.DR-10.0.somatic.maf -partition general -tumor_name CESC -NCI_skip 5 -cov CESC -Local_MAF mutationsTN_32_Cervical_cancer_no_cell_lines_just_SCC.maf -cores 12 -mem 5300 -time 24:00:00 -tumor_type_trinuc_sig CESC; cd ../
cd COAD; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf -partition general -tumor_name COAD -NCI_skip 5 -cov colon -Local_MAF mutationsTN_41_Metastatic_colon_cancer_to_liver__lung_and_elsewhere_JUST_primary.maf -cores 12 -mem 5300 -time 24:00:00 -tumor_type_trinuc_sig COAD; cd ../
cd ESCA; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.ESCA.mutect.7f8e1e7c-621c-4dfd-8fad-af07c739dbfc.DR-10.0.somatic.maf -partition general -tumor_name ESCA -NCI_skip 5 -cov ESCA -cores 12 -mem 5300 -time 24:00:00 -tumor_type_trinuc_sig Eso-AdenoCA; cd ../
cd GBM; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.GBM.mutect.da904cd3-79d7-4ae3-b6c0-e7127998b3e6.DR-10.0.somatic.maf -partition general -tumor_name GBM -NCI_skip 5 -cov GBM -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig GBM; cd ../
cd HNSC_HPVneg; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -cov HNSC -partition general -ready_MAF HNSC_HPV_neg.RData -tumor_name HNSC_HPVneg -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig HNSC; cd ../
cd HNSC_HPVpos; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -cov HNSC -partition general -ready_MAF HNSC_HPV_pos.RData -tumor_name HNSC_HPVpos -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig HNSC; cd ../
cd KIRC; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.KIRC.mutect.2a8f2c83-8b5e-4987-8dbf-01f7ee24dc26.DR-10.0.somatic.maf -partition general -tumor_name KIRC -NCI_skip 5 -cov KIRC -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig KIRC; cd ../
cd LGG; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.LGG.mutect.1e0694ca-fcde-41d3-9ae3-47cfaf527f25.DR-10.0.somatic.maf -partition general -tumor_name LGG -NCI_skip 5 -cov LGG -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig LGG; cd ../
cd LUAD; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf -Local_MAF mutationsTN_62_Lung_Adenocarcinoma_Yale___MD_Anderson.maf -partition general -tumor_name LUAD -NCI_skip 5 -cov lung -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig LUAD; cd ../
cd LUSC; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript  ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.LUSC.mutect.95258183-63ea-4c97-ae29-1bae9ed06334.DR-10.0.somatic.maf -Local_MAF mutationsTN_63_Lung_Squamous_Cell_Carcinoma_Yale___MD_Anderson.maf -NCI_skip 5 -partition general -tumor_name LUSC -cov lung -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig LUSC; cd ../
cd LIHC; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -cov LIHC -tumor_name LIHC -NCI_skip 5 -partition general -NCI_MAF TCGA.LIHC.mutect.a630f0a0-39b3-4aab-8181-89c1dde8d3e2.DR-10.0.somatic.maf -Local_MAF mutationsTN_29_Hepatocellular_Carcinoma.maf -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig LIHC; cd ../
cd OV; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.OV.mutect.b22b85eb-2ca8-4c9f-a1cd-b77caab999bd.DR-10.0.somatic.maf -partition general -tumor_name OV -NCI_skip 5 -cov OV -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig OV; cd ../
cd PAAD; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript  ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf -Local_MAF new_mutationsTN_26_Pancreatic_Cancer.maf -NCI_skip 5 -partition general -tumor_name PAAD -cov pancreas -cores 12 -mem 6000 -time 24:00:00 -tumor_type_trinuc_sig PAAD; cd ../
cd PRAD; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript  ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.PRAD.mutect.deca36be-bf05-441a-b2e4-394228f23fbe.DR-10.0.somatic.maf -NCI_skip 5 -partition general -tumor_name PRAD -cov PRAD -cores 12 -mem 6000  -time 24:00:00 -tumor_type_trinuc_sig PRAD; cd ../
cd READ; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript  ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf -NCI_skip 5 -partition general -tumor_name READ -cov rectum  -cores 12 -mem 6000  -time 24:00:00 -tumor_type_trinuc_sig COAD; cd ../
cd SKCM_primary; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -cov SKCM -partition general -ready_MAF SKCM_NCI_Yale_prim.RData -tumor_name SKCMP -cores 12 -mem 6000  -time 24:00:00 -tumor_type_trinuc_sig SKCM; cd ../
cd SKCM_metastasis; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript ../selection_and_pop_level_cluster_submitter.R -cov SKCM -partition general -ready_MAF SKCM_NCI_Yale_met.RData -tumor_name SKCMM -cores 12 -mem 6000  -time 24:00:00 -tumor_type_trinuc_sig SKCM; cd ../
cd STAD; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript  ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.STAD.mutect.c06465a3-50e7-46f7-b2dd-7bd654ca206b.DR-10.0.somatic.maf -NCI_skip 5 -partition general -tumor_name STAD -cov stomach -cores 12 -mem 6000  -time 24:00:00 -tumor_type_trinuc_sig Stomach-AdenoCA; cd ../
cd THCA; ~/../../pi/townsend/vlc24/R_source/R-3.6.1/bin/Rscript  ../selection_and_pop_level_cluster_submitter.R -NCI_MAF TCGA.THCA.mutect.13999735-2e70-439f-a6d9-45d831ba1a1a.DR-10.0.somatic.maf -NCI_skip 5 -partition general -tumor_name THCA -cov THCA -Local_MAF mutationsTN_36_Anaplastic_thyroid_carcinoma.maf -cores 12 -mem 6000  -time 24:00:00 -tumor_type_trinuc_sig THCA; cd ../