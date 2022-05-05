library(data.table)
data_path = commandArgs(t=T)[1]
rep = commandArgs(t=T)[2]
genome_folder = commandArgs(t=T)[3]

folder = sprintf("%s/replicate%s",data_path, rep)

QTL_pos = read.csv(sprintf("%s/QTLpos_rep%s.csv", folder, rep), header = T)[,-1]
QTL_names = as.character(unlist(QTL_pos))

gwasResults = fread(sprintf("%s/gcta_gwas_res_preselect_LDclump_rep%s.clumped",folder, rep))[,1:6]

preselected_QTL = intersect(QTL_names,gwasResults$SNP)
print("number of QTL have been selected after LD clumping")
length(preselected_QTL)

write.csv(preselected_QTL, sprintf("%s/preselectedQTL_LDclumped_rep%s.csv",folder, rep))

### check if they are top 
SNP_rank = 1:nrow(gwasResults)
print("rank of preselected QTL")
SNP_rank[gwasResults$SNP %in% preselected_QTL]

# extract snps for following MegaLMM analysis 
SNPs_selected_LDclumped = gwasResults$SNP
write.table(SNPs_selected_LDclumped,file = sprintf('%s/SNPs_selected_LDclumped_rep%s.txt',folder, rep),row.names=F,quote=F,col.names = F)

