## ST allFactor RMSE for SNP and QTL 

data_path = commandArgs(t=T)[1]
nTraits = as.numeric(commandArgs(t=T)[2]) #2 or 18
nFactors = as.numeric(commandArgs(t=T)[3]) # 2 or 6 or 9
nSNPs = as.numeric(commandArgs(t=T)[4])
folder = commandArgs(t=T)[5]
n_repeat = as.numeric(commandArgs(t=T)[6])

library(ROCR)
library(data.table)
library(ggplot2)

print(sprintf("nTraits = %d", nTraits))
print(sprintf("nFactors = %d", nFactors))
print(sprintf("nSNPs = %d", nSNPs))
print(sprintf("folder = %s", folder))
print(sprintf("n_replicate = %d", n_repeat))

ST_allFactor_markerEffect_RMSE = function(nFactors, nTraits, n_replicates, nSNPs, folder){
  
  Simulation_Folder = sprintf("%s/%s/%dfactor_%dtrait_%dsnp_0.95h2F", data_path, folder, nFactors, nTraits, nSNPs)
  print(Simulation_Folder)
  
  p = 2000
  
  RMSE_SNP = rep(0,n_replicates)
  RMSE_QTL = rep(0,n_replicates)
  
  for (r in 1:n_replicates){
    
    runID = sprintf("%s/MegaBayesC_%dFactors_%dTraits_rep%d", Simulation_Folder,nFactors,nTraits,r)
    
    mrk_eff_sample = fread(sprintf("%s/ST_rep%d/_marker_effects_yNA.txt", Simulation_Folder,r))
    mrk_eff_mean = colMeans(mrk_eff_sample)
    
    # true QTL position
    trueQTL_position = as.data.frame(fread(sprintf("%s/QTLpos_rep%d.csv",Simulation_Folder,r))[,-1])
    trueQTL_effect = as.data.frame(fread(sprintf("%s/QTLeffects_rep%d.csv",Simulation_Folder,r))[,-1])
    
    TrueLambda = as.matrix(fread(sprintf("%s/TrueLambda_%dfactors_%dtraits_rep%d.csv", Simulation_Folder, nFactors, nTraits, r))[,-1])
    
    true.ME.pos = unique(unlist(trueQTL_position))
    
    # true Marker effects
    B2F_true = matrix(0, p, nFactors)
    for(k in 1:nFactors){
      B2F_true[trueQTL_position[[k]],k] = trueQTL_effect[[k]]
    }
    B2F_star_true = B2F_true%*% TrueLambda[,1]
    
    # RMSE for QTL
    sse_qtl =sum((mrk_eff_mean[true.ME.pos]-B2F_star_true[true.ME.pos])^2)
    RMSE_QTL[r] = sqrt(sse_qtl/(length(true.ME.pos)))
    
    # RMSE for SNP
    sse_snp =sum(mrk_eff_mean[-true.ME.pos]^2)
    RMSE_SNP[r] = sqrt(sse_snp/(p-length(true.ME.pos)))
    
    print(sprintf("replicate %d is done.", r))
    
  }
  
  names(RMSE_QTL) = paste0("repeat",1:n_replicates)
  write.csv(RMSE_QTL, sprintf("%s/ST_RMSE_QTL.csv",Simulation_Folder))
  
  names(RMSE_SNP) = paste0("repeat",1:n_replicates)
  write.csv(RMSE_SNP, sprintf("%s/ST_RMSE_SNP.csv",Simulation_Folder))
  
  return(list(RMSE_SNP=RMSE_SNP, RMSE_QTL=RMSE_QTL))
}

ST_allFactor_markerEffect_RMSE(nFactors, nTraits, n_repeat, nSNPs, folder)

