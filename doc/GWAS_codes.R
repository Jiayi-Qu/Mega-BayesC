posterior_paths = "/Users/apple/Desktop/JiayiQu/UCD_PhD/SparseFactorModel/MegaLMM_data_analysis/Wheat/server_MCMC/MegaLMM_Krause_BayesC_scenario1_10/Posterior/"
samples_number = seq(50,4000,by=50)
predAccuracy = c()

for(s in samples_number){
  print(s)
  MegaLMM_state$current = readRDS(paste0(posterior_paths,"Lambda_", s,".rds"))
  MegaLMM_state$Posterior$Lambda =  readRDS(paste0(posterior_paths,"Lambda_", s,".rds"))
  MegaLMM_state$Posterior$B2_F = readRDS(paste0(posterior_paths,"B2_F_",s,".rds"))
  U = get_posterior_mean(MegaLMM_state,X2_F %*% B2_F %*% Lambda[,1])
  g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = U[nas,1]),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]
  print(g_cor)
  predAccuracy = c(predAccuracy, g_cor)
  rm(MegaLMM_state)
  MegaLMM_state = MegaLMM_state_base
  gc()
}

plot(seq(50,4000,by=50), predAccuracy, ylim = c(0,1), xlab = "iterations", type = "l")

library(ggplot2)
#### Window_based SNP = 1
#### Latent trait = 3
### Last iteration

X2 = MegaLMM_state1$data_matrices$X2_F
B2F = MegaLMM_state1$Posterior$B2_F[50,,]
Lambda = MegaLMM_state1$Posterior$Lambda[50,,]


Sample_WPPA_Single_SNP = function(X2, B2F, Lambda, latent_num, observed_trait_index){
  p = ncol(X2)
  ratio = rep(0,p)
  #### For SNP 1 (local variance for observed-trait 1 given 3 latent traits)
  for(snp_index in 1:p){
    local_var = as.matrix(X2[,snp_index]) %*% B2F[snp_index,1:latent_num] %*% Lambda[1:latent_num,observed_trait_index]
    latent_trait_global_var = X2 %*% B2F[,1:latent_num] %*% Lambda[1:latent_num,observed_trait_index]
    #all_global_var = X2 %*% B2F %*% Lambda[,observed_trait_index]
    snp_wppa = var(local_var)/var(latent_trait_global_var)
    ratio[snp_index] = snp_wppa
  }
  
  return(list(ratio = ratio))
}


p = ncol(X2)

WPPA_results1 = data.frame(matrix(nrow = 50, ncol = (p+1)))
names(WPPA_results1) = c("iteration", sapply(1:p,function(x) paste0("SNP",x)))

for(i in 1:50){
  B2F = MegaLMM_state1$Posterior$B2_F[i,,]
  Lambda = MegaLMM_state1$Posterior$Lambda[i,,]
  res = Sample_WPPA_Single_SNP(X2, B2F, Lambda, latent_num=1, observed_trait_index=1)
  WPPA_results1[i, ] = c(i, res$ratio)
  print(sprintf("iteration %d finished.", i))
}

WPPA_results1 = WPPA_results1[,-1]
wppa_threshold = 1/p
wppa_identity1 = as.numeric(WPPA_results1 > wppa_threshold)


plot_df = data.frame(SNP = 1:p, WPPA = colMeans(wppa_identity1))

highlight_df <- plot_df %>% filter(WPPA >0.65)
max_df <- highlight_df[which.max(highlight_df$WPPA),]

png(paste0("WPPA_10factor.jpg"), width = 7, height = 7, units = 'in', res = 300)
ggplot(data = plot_df, aes(x=SNP,y=WPPA)) + 
  geom_point(alpha=0.3) + ylim(0:1) + 
  geom_point(data=highlight_df, aes(x=SNP,y=WPPA), color='red',size=2) +
  geom_label(
    label=rownames(max_df), 
    x=max_df$SNP,
    y=(max_df$WPPA+0.1),
    label.padding = unit(0.55, "lines"), # Rectangle size around label
    label.size = 0.35,
    color = "black",
    fill="#69b3a2"
  ) + ggtitle("Manhattan Plot ( factor)")
dev.off()







