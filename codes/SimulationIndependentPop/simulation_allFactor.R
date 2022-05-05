library(data.table)
##  Here, all factors are assumed to have effects on the focal trait

#simulation_path = "/Users/apple/Desktop/JiayiQu/UCD_PhD/SparseFactorModel/MegaLMM_data_analysis/Wheat"
simulation_path = commandArgs(t=T)[1]
#code_path = "/Users/apple/Desktop/JiayiQu/UCD_PhD/SparseFactorModel/MegaLMM_data_analysis/Wheat"
code_path = commandArgs(t=T)[2]

#MeanLambdaValue = as.numeric(commandArgs(t=T)[4])

nTraits = as.numeric(commandArgs(t=T)[3]) #2 or 18
nFactors = as.numeric(commandArgs(t=T)[4]) # 2 or 6 or 9
#nTraits = 18
#nFactors = 9

TrueLambda = matrix(0, nrow = nFactors, ncol = nTraits)

nInd = 5000

#numberofsnp = 10
numberofsnp = as.numeric(commandArgs(t=T)[5])
#replicate = 1
replicate = as.numeric(commandArgs(t=T)[6])

set.seed(replicate)


NonZerosPerRow = nTraits/nFactors

#Lambda_values = runif(nTraits,MeanLambdaValue-0.1,MeanLambdaValue+0.1)
Lambda_values = rep(1,nTraits)
for (r in 1:nFactors){
  TrueLambda[r, (NonZerosPerRow*(r-1)+1):(NonZerosPerRow*r)] = Lambda_values[(NonZerosPerRow*(r-1)+1):(NonZerosPerRow*r)]
}

# all factors are linked to focal trait
TrueLambda[2:nFactors, 1] = 1


colnames(TrueLambda) = c("focal",paste0("secondary",1:(nTraits-1)))
rownames(TrueLambda) = paste0("factor",1:nFactors)

write.csv(TrueLambda, sprintf("%s/TrueLambda_%dfactors_%dtraits_rep%d.csv",simulation_path,nFactors,nTraits,replicate))


nMarkers = 2000

h2_R = 0.9
h2_F = 0.95
#h2_R = as.numeric(commandArgs(t=T)[7])
#h2_R = as.numeric(commandArgs(t=T)[8])


set.seed(replicate)
X2 = matrix(sample(c(0,1,2), nInd*nMarkers, replace = T), nrow = nInd, ncol = nMarkers)
colnames(X2) = paste0("SNP", 1:nMarkers)
rownames(X2) = paste0("Ind", 1:nInd)
write.csv(X2, sprintf("%s/X2_rep%d.csv",simulation_path,replicate))





### Model = (X_2 B_2F + E_F) Lambda + E_R
### For each observed trait, we generate the corresponding phenotypic values


K = nFactors 

effective_factor = 1:K
Y = matrix(0, nrow = nInd, ncol=nTraits)

nQTL = rep(numberofsnp,K)

### assume each latent trait's heritability = var(X_2 b_2F)/(var(X_2 b_2F) + var(e_f))
### assume each latent trait's heritability = 0.95
h2F = rep(h2_F, K)

### Sample the QTL ID and QTL effects

#QTL_effects = lapply(nQTL,  function(x) rnorm(x))

QTL_effects = lapply(nQTL, function(x) rep(0.1,x))

names(QTL_effects) = paste0("factor", effective_factor)
write.csv(as.data.frame(QTL_effects), sprintf("%s/QTLeffects_rep%d.csv",simulation_path,replicate))


QTL_pos = lapply(nQTL,  function(x) sample(1:nMarkers, x, replace = FALSE))
names(QTL_pos) = paste0("factor", effective_factor)
write.csv(as.data.frame(QTL_pos), sprintf("%s/QTLpos_rep%d.csv",simulation_path,replicate))
print("QTL information is saved.")


### Simulate Breeding value (nIndxn_factors matrix)
BVs = sapply(1:K, function(x) X2[, QTL_pos[[x]]] %*% QTL_effects[[x]])
colnames(BVs) = paste0("factor", effective_factor)
rownames(BVs) = rownames(X2)

# variance explained per factor before normalization
var_explained_per_factor = apply(BVs,2,var)
write.csv(var_explained_per_factor, sprintf("%s/var_explained_per_factor_rep%d.csv",simulation_path,replicate))


# normalize BVs
#BVs = apply(BVs, 2, function(x) x/sd(x))


### Simulate E_F
vare_f = (1-h2F)/h2F * var_explained_per_factor
EF = sapply(1:K, function(x) rnorm(nInd) * sqrt(vare_f[x]))

### Simulate F 
F_matrix = BVs +  EF

print(sprintf("The simulated h2F for factor %d is %f",1:K, 1/(apply(EF,2,var)+1)))



### Simulate er (h2_er = var(Fxlambda_i)/var(Fxlambda_i + er) = 0.9)
FXLambda = F_matrix %*% as.matrix(TrueLambda)
FXLambda_var = apply(FXLambda,2,var)
vare_r = (1-h2_R)/h2_R * FXLambda_var
er = sapply(1:nTraits, function(x) rnorm(nInd) * sqrt(vare_r[x]))

### Generates the  observed_traits in Y 
Y = FXLambda + er
colnames(Y) = c("focal",paste0("secondary",1:(nTraits-1)))

Y_var = apply(Y,2,var)
Bstar = BVs %*% as.matrix(TrueLambda)
Bstar_var = apply(Bstar,2,var)
h2_final = Bstar_var/Y_var
h2_final

write.csv(h2_final, sprintf("%s/h2_final_%dFactors_%dTraits_rep%d.csv", simulation_path,nFactors, nTraits,replicate))
print("The h2_final is")
print(h2_final)

write.csv(Y, sprintf("%s/Y_%dFactors_%dTraits_rep%d.csv",simulation_path,nFactors, nTraits,replicate))
print("Y matrix is saved.")

source(sprintf("%s/data_prep_Krause.R",code_path))

source(sprintf("%s/MegaLMM_Krause_geno_BayesC.R",code_path))


