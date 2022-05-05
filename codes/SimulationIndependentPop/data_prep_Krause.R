library(data.table)
library(rrBLUP)
library(tidyr)
library(lme4)
library(MCMCglmm)
library(lme4qtl)
library(MegaLMM)

set.seed(replicate)

BLUEs = fread(sprintf("%s/Y_%dFactors_%dTraits_rep%d.csv",simulation_path,nFactors,nTraits,replicate),data.table=F)
colnames(BLUEs)[1] ="GID"
BLUEs$GID = as.character(BLUEs$GID)

data = data.frame(GID = BLUEs$GID, BLUE = BLUEs$focal)

print(replicate)

