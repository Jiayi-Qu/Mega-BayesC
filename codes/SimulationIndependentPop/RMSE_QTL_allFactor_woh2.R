data_path = commandArgs(t=T)[1]
nTraits = as.numeric(commandArgs(t=T)[2]) #2 or 18
nFactors = as.numeric(commandArgs(t=T)[3]) # 2 or 6 or 9
nSNPs = as.numeric(commandArgs(t=T)[4])
folder = commandArgs(t=T)[5]
n_repeat= as.numeric(commandArgs(t=T)[6])

library(ROCR)
library(data.table)
library(ggplot2)

print(sprintf("nTraits = %d", nTraits))
print(sprintf("nFactors = %d", nFactors))
print(sprintf("nSNPs = %d", nSNPs))
print(sprintf("folder = %s", folder))
print(sprintf("n_replicate = %d", n_repeat))

allFactor_wo_h2_QTLEffect_RMSE = function(nFactors, nTraits, n_replicates, nSNPs, folder = folder){
  
  Simulation_Folder = sprintf("%s/%s/%dfactor_%dtrait_%dsnp_0.95h2F", data_path, folder, nFactors, nTraits, nSNPs)
  print(Simulation_Folder)
  
  p = 2000
  K = 10
  
  RMSE_QTL = rep(0,n_replicates)
  
  for (r in 1:n_replicates){
    
    runID = sprintf("%s/MegaBayesC_%dFactors_%dTraits_rep%d", Simulation_Folder,nFactors,nTraits,r)
    
    # posterior mean of marker effects
    B2F_PosteriorMean = read.csv(sprintf('%s/Posterior/Posterior_B2_F.csv',runID))
    B2F_PosteriorMean = B2F_PosteriorMean[,-1]
    
    
    # true QTL position
    trueQTL_position = as.data.frame(fread(sprintf("%s/QTLpos_rep%d.csv",Simulation_Folder,r))[,-1])
    trueQTL_effect = as.data.frame(fread(sprintf("%s/QTLeffects_rep%d.csv",Simulation_Folder,r))[,-1])
    
    
    # estimated and true marker effects 
    estLambda = as.matrix(fread(sprintf('%s/Posterior/Posterior_Lambda.csv',runID))[,-1])
    TrueLambda = as.matrix(fread(sprintf("%s/TrueLambda_%dfactors_%dtraits_rep%d.csv", Simulation_Folder, nFactors, nTraits, r))[,-1])
    
    B2F_star = as.matrix(B2F_PosteriorMean) %*% estLambda[,1]   ### focal trait index = 1 
    
    true.ME.pos = trueQTL_position$factor1
    
    # true Marker effects
    B2F_true = matrix(0, p, nFactors)
    k = 1 
    B2F_true[trueQTL_position[[k]],k] = trueQTL_effect[[k]]
    B2F_star_true = B2F_true%*% TrueLambda[,1]
    
    
    # RMSE for QTL
    sse_qtl =sum((B2F_star[true.ME.pos]-B2F_star_true[true.ME.pos])^2)
    RMSE_QTL[r] = sqrt(sse_qtl/(length(true.ME.pos)))
    
    
    print(sprintf("replicate %d is done.", r))
    
  }
  
  names(RMSE_QTL) = paste0("repeat",1:n_replicates)
  write.csv(RMSE_QTL, sprintf("%s/RMSE_QTL.csv",Simulation_Folder))
  
  return(list(RMSE_QTL=RMSE_QTL))
}

allFactor_wo_h2_QTLEffect_RMSE(nFactors, nTraits, n_repeat, nSNPs, folder = folder)
