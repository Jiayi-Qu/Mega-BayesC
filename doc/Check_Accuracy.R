path = "/Users/apple/Desktop/JiayiQu/UCD_PhD/SparseFactorModel/MegaLMM_data_analysis/Wheat/Results_1415_OF_Bgcor/server_results/"
files = list.files(path)
results = c()
for(f in files){
  d = fread(paste0(path,f))
  results = rbind(results, d)
}

results = results[,-1]
mean(results[results$Method == "MegaLMM_BayesC_s1",]$g_cor)
mean(results[results$Method == "MegaLMM_BayesC_s4",]$g_cor)
mean(results[results$Method == "MegaLMM_BayesC_s5",]$g_cor)
mean(results[results$Method == "MegaLMM_BayesC_s6",]$g_cor)







