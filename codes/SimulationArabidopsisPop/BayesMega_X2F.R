library(data.table)
library(stringr)
library(plyr)  
data_path = commandArgs(t=T)[1]
rep = commandArgs(t=T)[2]
genome_folder = commandArgs(t=T)[3]

folder = sprintf("%s/replicate%s",data_path, rep)

preselected_matrix = read.table(sprintf("%s/BayesMega_preselected_Matrix_LDclumped_rep%s.raw", folder, rep), header = T)
rownames(preselected_matrix) = preselected_matrix$IID
preselected_matrix = preselected_matrix[, -(1:6)]

#Imputation (mean)
for (j in 1:ncol(preselected_matrix)) {
  preselected_matrix[, j] = ifelse(is.na(preselected_matrix[, j]), mean(preselected_matrix[, j], na.rm = TRUE), preselected_matrix[, j])
}

SNP_names = colnames(preselected_matrix)
SNP_names = str_split_fixed(SNP_names, "_",2)[,1]
colnames(preselected_matrix) = SNP_names

write.csv(preselected_matrix,sprintf("%s/geno_file_LDclumped_rep%s.csv",folder, rep))

M_qtl = read.csv(sprintf("%s/QTL_genotype_imputed_GeneExpInd_rep%s.csv", folder, rep))
rownames(M_qtl) = M_qtl[,1]
M_qtl = M_qtl[,-1]


selected_QTL = intersect(colnames(preselected_matrix), colnames(M_qtl))
non_selected_QTL = colnames(M_qtl)[!colnames(M_qtl) %in% selected_QTL]
M_reducedQTL = preselected_matrix[,!colnames(preselected_matrix) %in% selected_QTL]

CorbtwQTLASNP = matrix(0, nrow = length(non_selected_QTL), ncol = 3)
colnames(CorbtwQTLASNP) = c("QTL","max_cor","SNP")

i = 1 
for (qtl in non_selected_QTL){
  print(sprintf("QTL %s.", qtl))
  CorbtwQTLASNP[i, "QTL"] = qtl
  r2 = cor(M_qtl[,qtl], M_reducedQTL)^2
  r2 = t(r2)
  max_index = which.max(r2[,1])
  max_value = max(r2[,1])
  CorbtwQTLASNP[i, "max_cor"] = max_value
  max_snp_name = names(r2[max_index,])
  CorbtwQTLASNP[i, "SNP"] = max_snp_name
  print(max_snp_name)
  print(max_value)
  i = i+1
}

write.csv(CorbtwQTLASNP, sprintf("%s/CorbtwQTLASNP_LDclumped_rep%s.csv",folder, rep))


