library(stringr)
library(data.table)
library(MegaLMM)
library(plyr)

data_path = commandArgs(t=T)[1]
rep = commandArgs(t=T)[2]
genome_folder = commandArgs(t=T)[3]

folder = sprintf("%s/replicate%s",data_path, rep)

# simulated Lambda
nFactor =  10
nTrait = nFactor*20 + 1

True_sim_Lambda  = matrix(0, nrow = nFactor, ncol = 20*nFactor)
for (i in 1:nFactor){
  True_sim_Lambda[i, ((i-1)*20+1):(i*20)] = 1
}
True_sim_Lambda = cbind(rep(0.5, nFactor), True_sim_Lambda)
True_sim_Lambda[1:10,11:20]
write.csv(True_sim_Lambda, sprintf("%s/True_sim_Lambda_rep%s.csv",folder, rep))

# Scenario 1.1 non-pleiotropy, different QTL affect different factors (no linked traits)
F_h2 = rep(1, 10) # assume E_F = 0
resid_h2 = c(0.6, rep(0.8,20*nFactor))

# Extract 20 QTL from the true QTL matrix
Raw_df = read.table(sprintf("%s/QTL_sim_info_rep%s.raw", folder, rep), header = T)
rownames(Raw_df) = Raw_df$IID
Raw_df = Raw_df[, -(1:6)]
Raw_df[1:10,1:10]

for (j in 1:ncol(Raw_df)) {
  Raw_df[, j] = ifelse(is.na(Raw_df[, j]), mean(Raw_df[, j], na.rm = TRUE), Raw_df[, j])
}


W_qtl = Raw_df[!duplicated(as.list(Raw_df))]
print("dim of M_qtl after rm duplicated columns (should be 1003 x 20)")
dim(W_qtl) 

set.seed(as.numeric(rep))

SNP_names = colnames(W_qtl)
SNP_names = str_split_fixed(SNP_names, "_",2)[,1]
colnames(W_qtl) = SNP_names
print("sd of QTL genotypic values:")
apply(W_qtl, 2, sd)

QTL_names = sample(SNP_names)

write.csv(W_qtl, sprintf("%s/QTL_genotype_imputed_rep%s.csv",folder, rep))


# Simulation process
K = nFactor 
nInd = nrow(W_qtl)

numberofqtl = 2
effective_factor = 1:K
Y = matrix(0, nrow = nInd, ncol=nTrait)
nQTL = rep(numberofqtl,K)

# The h2 of latent trait var(X2B2F)/var(X2B2F+EF)
h2F = as.numeric(F_h2)
names(h2F) = paste0("factor",effective_factor)

QTL_effects = lapply(nQTL, function(x){
  positive = sample(c(T,F), 1)
  if (positive) {
    runif(x, min = 3, max = 5)
  }else{
    runif(x, min = -5, max = -3)
  }
})
head(QTL_effects)

##### generate the information for  QTL 
names(QTL_effects) = paste0("factor", effective_factor)
write.csv(as.data.frame(QTL_effects), sprintf("%s/QTLeffects_rep%s.csv",folder,rep))


# every (i-1)*numberofqtl:i*numberofqtl is the QTL for factior i 
QTL_pos = lapply(1:K,  function(x){
  QTL_names[((x-1)*numberofqtl + 1):(x*numberofqtl)]
})
names(QTL_pos) = paste0("factor", effective_factor)
head(QTL_pos)
write.csv(as.data.frame(QTL_pos), sprintf("%s/QTLpos_rep%s.csv",folder, rep))
print("QTL information is saved.")


### Simulate Breeding value (nIndxn_factors matrix)
BVs = sapply(1:K, function(x) as.matrix(W_qtl[, QTL_pos[[x]]]) %*% QTL_effects[[x]])
colnames(BVs) = paste0("factor", effective_factor)
rownames(BVs) = rownames(W_qtl)
BVs[1:10,1:10]

# variance explained per factor before normalization
var_explained_per_factor = apply(BVs,2,var)
write.csv(var_explained_per_factor, sprintf("%s/var_explained_per_factor_rep%s.csv",folder, rep))
var_explained_per_factor

### Simulate E_F
vare_f = (1-h2F)/h2F * var_explained_per_factor
EF = sapply(1:K, function(x) rnorm(nInd) * sqrt(vare_f[x]))

### Simulate F
F_matrix = BVs +  EF
simulated_h2F = apply(BVs,2,var)/apply(F_matrix,2,var)
write.csv(simulated_h2F, sprintf("%s/simulated_h2F_rep%s.csv",folder, rep))

# normalize F by the equation given by Runcie (corrected)
EF_norm = sapply(1:K, function(x) rnorm(nInd) * sqrt(1-simulated_h2F[x]))
BV_norm = sapply(1:K, function(x) BVs[,x]/sqrt(var_explained_per_factor[x]/simulated_h2F[x]))
F_norm = EF_norm + BV_norm
apply(F_norm, 2, var)

FXLambda = F_norm %*% as.matrix(True_sim_Lambda)
var_per_trait = apply(FXLambda,2,var)

### Simulate Er 
# (1-resid_h2)*sigma^2_F/resid_h2, resid_h2 = sigma^2_F/(sigma^2_F+sigma^2_R)
residual_sd = sqrt((1-resid_h2)/resid_h2 * var_per_trait)
head(residual_sd)
ER =  sapply(1:nTrait, function(x) rnorm(nInd, sd = residual_sd[x]))

### Generates the observed_traits in Y

Y = FXLambda + ER
colnames(Y) = c("focal", paste0("secondary",1:(20*nFactor)))

# update the marker effects 
QTL_effect_updated = lapply(1:K,  function(x){
  QTL_effects[[x]]/sqrt(var_explained_per_factor[x]/simulated_h2F[x])
})
names(QTL_effect_updated) = paste0("factor", effective_factor)
write.csv(as.data.frame(QTL_effect_updated), sprintf("%s/QTLeffects_updated_rep%s.csv",folder, rep))


lambda.1 = True_sim_Lambda[,1]
# true marker effects
QTL_effect_updated = read.csv(sprintf("%s/QTLeffects_updated_rep%s.csv",folder, rep))[,-1]
QTL_pos = as.data.frame(QTL_pos)
QTLs = as.character(unlist(QTL_pos))
QTLs
# calculate the true focal QTL effects 
true_focal_QTL_effects_list = lapply(QTLs, function(x){
  factor_boolean = apply(QTL_pos ==  x, 2, sum)
  factor_name = names(factor_boolean[factor_boolean == 1])
  factor_index = as.numeric(str_replace(factor_name, "factor",""))
  true_focal_effect = QTL_effect_updated[QTL_pos ==  x] * lambda.1[factor_index]
  true_focal_effect 
})
names(true_focal_QTL_effects_list) = QTLs
true_focal_QTL_effects = unlist(true_focal_QTL_effects_list)
write.csv(true_focal_QTL_effects, sprintf("%s/true_focal_QTL_effects_rep%s.csv", folder, rep))


#plot(unlist(QTL_effects), true_focal_QTL_effects)
# look at the distribution of effect sizes for focal trait
#hist(true_focal_QTL_effects)
num_qtl = 20
### variance explained by each QTL for focal trait
var_qtl_explained = sapply(1:num_qtl, function(x) {
  marker = true_focal_QTL_effects[x]
  markerID = names(marker)
  X2 = W_qtl[,markerID]
  local_EBV = X2 * marker
  var = var(local_EBV)
  names(var) = markerID
  var
})

var_per_QTL = var_qtl_explained/var(Y[,1]) 
write.csv(var_per_QTL, sprintf("%s/var_per_QTL_rep%s.csv", folder, rep))

print("number of QTL explained more than 1% variance: (20 total)")
sum(var_per_QTL * 100 > 0.1)

png(sprintf("%s/hist_var_per_QTL_rep%s.jpg",folder, rep), width = 7, height = 7, units = 'in', res = 300)
print(hist(var_per_QTL*100))
dev.off()

Y_var = apply(Y,2,var)
XBstar = BV_norm %*% as.matrix(True_sim_Lambda)
XBstar_var = apply(XBstar,2,var)
h2_final = XBstar_var/Y_var
head(h2_final)
#hist(h2_final)

write.csv(h2_final, sprintf("%s/h2_final_rep%s.csv", folder, rep))

write.csv(Y, sprintf("%s/Y_rep%s.csv",folder, rep))

print("Y matrix is saved.")

## generate FT10_sim file
FT10_sim = data.frame(accession_id = rownames(Y), phenotype_value = Y[,1])
write.csv(FT10_sim, sprintf("%s/FT10_sim_all_ind_rep%s.csv",folder, rep), row.names = F, quote = F)
head(FT10_sim)

## only select individuals having GEXP values 
GEXP_ind = read.table(sprintf("%s/gene_exp_ind_data.txt",genome_folder))
subW_gexp = W_qtl[as.character(GEXP_ind$V1),]
dim(subW_gexp)
write.csv(subW_gexp, sprintf("%s/QTL_genotype_imputed_GeneExpInd_rep%s.csv",folder, rep))


FT10_sim_GEXP = FT10_sim[match(GEXP_ind$V1,FT10_sim$accession_id),]
write.csv(FT10_sim_GEXP, sprintf("%s/FT10_sim_gexp_ind_rep%s.csv",folder, rep), row.names = F,quote=F)
head(FT10_sim_GEXP)

Y = as.data.frame(Y)
GEXP_trait = Y[match(GEXP_ind$V1,rownames(Y)),-1]
write.csv(GEXP_trait, sprintf("%s/GEXP_trait_rep%s.csv",folder, rep))

FT10_sim_preselect = FT10_sim[-match(GEXP_ind$V1,FT10_sim$accession_id),]
FT10_sim_preselect$accession_id = as.numeric(as.character(FT10_sim_preselect$accession_id))
head(FT10_sim_preselect)
pheno_gcta_sim = data.frame(FID=FT10_sim_preselect$accession_id,IID = FT10_sim_preselect$accession_id,pheno = FT10_sim_preselect$phenotype_value)
write.table(pheno_gcta_sim, sprintf("%s/pheno_gcta_sim_preselect.phen",folder), row.names = F, col.names = F)

# phen file for analysis using all ind
pheno_gcta_sim_all_ind = data.frame(FID=FT10_sim$accession_id,IID = FT10_sim$accession_id,pheno = FT10_sim$phenotype_value)
write.table(pheno_gcta_sim_all_ind, sprintf("%s/pheno_gcta_sim_all_ind.phen",folder), row.names = F, quote = F, col.names = F)








