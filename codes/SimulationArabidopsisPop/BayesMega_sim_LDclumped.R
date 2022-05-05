library(data.table)
library(rrBLUP)
library(tidyr)
library(lme4qtl)
library(lme4)
library(MegaLMM)

set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()-1)

foldid = 123
if(is.na(foldid)) foldid = 1
set.seed(foldid)

analysis_path = commandArgs(t=T)[1]
data_path = commandArgs(t=T)[2]
pop = as.numeric(commandArgs(t=T)[3]

runID = sprintf('%s/replicate%d/BayesMega_rep%d', analysis_path, pop, pop)

run_parameters = MegaLMM_control(
  which_sampler = list(Y=4,F=4),
  run_sampler_times = 1,
  drop0_tol = 1e-10,
  scale_Y = F,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
  h2_divisions = 5, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
  h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
  burn = 10000,  # number of burn in samples before saving posterior samples
  K = 30, # number of factors
  thin = 100
)


pheno = sprintf("%s/replicate%d/FT10_sim_gexp_ind_rep%d.csv",data_path, pop, pop)
pheno = read.csv(pheno)

data = data.frame(FT = pheno$phenotype_value,
                  ID = as.character(pheno$accession_id))

#GeneExp = fread(sprintf("%s/At_GE_VST.csv",data_path))
GeneExp = fread(sprintf("%s/replicate%d/GEXP_trait_rep%d.csv",data_path, pop, pop))
GeneExp$V1 = as.character(GeneExp$V1)
GeneExp = GeneExp[match(data$ID, GeneExp$V1),-1]
Y = cbind(data$FT, GeneExp)
rownames(Y) = data$ID
dim(Y)


priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5,   nu = 10),      # Prior variance of trait residuals after accounting for fixed effects and factors
  tot_F_var = list(V = 18/20, nu = 10000),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
 Lambda_prior = list(
    sampler = sample_Lambda_prec_BayesC, # function that implements the horseshoe-based Lambda prior described in Runcie et al 2020. See code to see requirements for this function.
    Lambda_df = 4,
    Lambda_scale = 1,
    delta_1 = list(shape = 20, rate = 1/2),
    delta_2 = list(shape = 3, rate = 1),
    fixed_pi = NULL,
    delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
  ),
  B2_prior = list(
    sampler = sample_B2_prec_BayesC,
    fixed_pi = NULL,
    B2_F_df = 4,
    B2_F_scale = 1,
    B2_R_df = 4,
    B2_R_scale = 1
  ),
  h2_priors_resids_fun = function(h2s,n) 1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
)

# use genotype of preselected SNPs in X
X = fread(sprintf("%s/replicate%d/geno_file_LDclumped_rep%d.csv",data_path, pop, pop),data.table=F, h=T)
rownames(X) = X[,1]
X = as.matrix(X[,-1])
dim(X)


MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                                    ~ 1 + (1|ID),
                                    extra_regressions = list(X = X,resids = F,factors=T),
                                    data=data,         # the data.frame with information for constructing the model matrices
                                    run_parameters=run_parameters,
                                    run_ID = runID
)

maps = make_Missing_data_map(MegaLMM_state,2,verbose=T)
MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)
MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state_base = MegaLMM_state
MegaLMM_state = initialize_MegaLMM(MegaLMM_state)
MegaLMM_state = clear_Posterior(MegaLMM_state)

MegaLMM_state$Posterior$posteriorSample_params = c("Lambda","F", "B2_F", "Lambda_pi", "B2_F_pi", "F_h2","tot_Eta_prec", "resid_h2")
MegaLMM_state = clear_Posterior(MegaLMM_state)

n_iter =500
for(i  in 1:100) {
  print(sprintf('Run %d',i))
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)  
 
  MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
  print(MegaLMM_state) # print status of current chain
  plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
  # set of commands to run during burn-in period to help chain converge
  if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i <= 20) {
    MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    print(MegaLMM_state$run_parameters$burn)
  }
}

X2_used = MegaLMM_state$data_matrices$X2_F
write.csv(X2_used, sprintf('%s/X2F.csv',runID))
dim(X2_used)


MegaLMM_state = MegaLMM_state_base
MegaLMM_state$Posterior = readRDS(sprintf('%s/Posterior/Posterior_base.rds',runID))

MegaLMM_state$Posterior$Lambda = load_posterior_param(MegaLMM_state,'Lambda')
MegaLMM_state$Posterior$B2_F = load_posterior_param(MegaLMM_state,'B2_F')
MegaLMM_state$Posterior$Lambda_pi = load_posterior_param(MegaLMM_state,'Lambda_pi')
MegaLMM_state$Posterior$B2_F_pi = load_posterior_param(MegaLMM_state,'B2_F_pi')
MegaLMM_state$Posterior$F_h2 = load_posterior_param(MegaLMM_state,'F_h2')
MegaLMM_state$Posterior$tot_Eta_prec = load_posterior_param(MegaLMM_state,'tot_Eta_prec')
MegaLMM_state$Posterior$resid_h2 = load_posterior_param(MegaLMM_state,'resid_h2')

saveRDS(MegaLMM_state$Posterior$Lambda,sprintf('%s/Posterior/Posterior_Lambda_samples.rds',runID))
saveRDS(MegaLMM_state$Posterior$B2_F,sprintf('%s/Posterior/Posterior_B2F_samples.rds',runID))
saveRDS(MegaLMM_state$Posterior$Lambda_pi,sprintf('%s/Posterior/Posterior_Lambda_pi_samples.rds',runID))
saveRDS(MegaLMM_state$Posterior$B2_F_pi,sprintf('%s/Posterior/Posterior_B2_F_pi_samples.rds',runID))
saveRDS(MegaLMM_state$Posterior$F_h2,sprintf('%s/Posterior/Posterior_F_h2_samples.rds',runID))
saveRDS(MegaLMM_state$Posterior$tot_Eta_prec,sprintf('%s/Posterior/Posterior_tot_Y_prec_samples.rds',runID))
saveRDS(MegaLMM_state$Posterior$resid_h2,sprintf('%s/Posterior/Posterior_resid_h2_samples.rds',runID))

print("RDS SAVED")

write.csv(get_posterior_mean(MegaLMM_state, Lambda), sprintf('%s/Posterior/Posterior_Lambda.csv',runID))
write.csv(get_posterior_mean(MegaLMM_state, B2_F), sprintf('%s/Posterior/Posterior_B2_F.csv',runID))
write.csv(get_posterior_mean(MegaLMM_state, Lambda_pi), sprintf('%s/Posterior/Posterior_Lambda_pi.csv',runID))
write.csv(get_posterior_mean(MegaLMM_state, B2_F_pi), sprintf('%s/Posterior/Posterior_B2_F_pi.csv',runID))
write.csv(get_posterior_mean(MegaLMM_state, F_h2), sprintf('%s/Posterior/Posterior_F_h2.csv',runID))
write.csv(get_posterior_mean(MegaLMM_state, tot_Eta_prec), sprintf('%s/Posterior/Posterior_tot_Y_prec.csv',runID))
write.csv(get_posterior_mean(MegaLMM_state, resid_h2), sprintf('%s/Posterior/Posterior_resid_h2.csv',runID))




